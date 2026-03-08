# eFISHent Hosted Service — MVP Plan

A minimal hosted version of eFISHent where users submit probe design jobs through a web interface and receive results via email/download.

## Architecture

```
┌─────────────────────────────────────────────────────────────┐
│                        Users                                │
│              (submit gene name + parameters)                │
└──────────────────────┬──────────────────────────────────────┘
                       │
                       ▼
┌──────────────────────────────────────────────────────────────┐
│  Always-On Frontend (t4g.small)                              │
│  ┌────────────────────────────────────────────────────────┐  │
│  │  FastAPI app                                          │  │
│  │  • Web form for job submission                        │  │
│  │  • Job status page / polling                          │  │
│  │  • Results download                                   │  │
│  │  • SQLite or DynamoDB for job metadata                │  │
│  └──────────────────────┬─────────────────────────────────┘  │
│                         │ submit job definition              │
│                         ▼                                    │
│  ┌────────────────────────────────────────────────────────┐  │
│  │  AWS Batch (submit API call)                          │  │
│  └──────────────────────┬─────────────────────────────────┘  │
└──────────────────────────────────────────────────────────────┘
                          │
                          ▼
┌──────────────────────────────────────────────────────────────┐
│  AWS Batch Compute Environment                               │
│  (Spot instances, c6i.2xlarge or c6i.4xlarge)                │
│                                                              │
│  ┌────────────────────────────────────────────────────────┐  │
│  │  Docker container (eFISHent + all dependencies)       │  │
│  │  • bowtie2, jellyfish, blastn, dustmasker, glpsol     │  │
│  │  • Fold binary (Linux)                                │  │
│  │  • Python + eFISHent package                          │  │
│  └────────────────────────────────────────────────────────┘  │
│                                                              │
│  Mounts:                                                     │
│  • S3: pre-built genome indices (read-only)                  │
│  • S3: job output bucket (write)                             │
└──────────────────────────────────────────────────────────────┘
```

### How It Works

1. User fills out a web form: gene name (or uploads FASTA), organism (human/mouse), preset (smFISH/MERFISH/DNA-FISH), optional advanced parameters
2. Frontend validates input, creates a job record, submits to AWS Batch
3. AWS Batch provisions a Spot instance, pulls the Docker image from ECR, mounts genome indices from S3
4. eFISHent runs inside the container (8–16 threads, 15–60 min typical)
5. Output files (`.fasta`, `.csv`, `.txt`, `.pdf`) are uploaded to S3
6. Frontend polls Batch job status; when complete, shows download link and/or sends email notification
7. Spot instance terminates automatically — no idle compute cost

### Why This Architecture

- **No idle compute cost** — Batch instances only exist while jobs run. A job that takes 30 minutes on a c6i.2xlarge spot costs ~$0.07
- **Simple frontend** — t4g.small (2 vCPU, 2 GB) is more than enough for a FastAPI app serving a handful of concurrent users
- **No orchestration complexity** — AWS Batch handles provisioning, queuing, retries, and cleanup. No need for Celery/Redis/Kubernetes
- **Spot interruption handling** — eFISHent jobs are idempotent (same input → same output). If a Spot instance is reclaimed, Batch retries the job automatically

---

## Pre-Built Resources

### Genome Indices (stored on S3, mounted at runtime)

| Organism | Genome | Bowtie2 Index | Jellyfish | GTF Parquet | Total |
|----------|--------|--------------|-----------|-------------|-------|
| Human | GRCh38 | ~3.5 GB | ~500 MB | ~200 MB | ~4.2 GB |
| Mouse | GRCm39 | ~3.0 GB | ~400 MB | ~150 MB | ~3.6 GB |
| **Total** | | | | | **~8 GB** |

Build these once on a larger instance, upload to S3. Jobs download/mount at startup (~30–60s on a fast instance).

### Docker Image

Single image containing:
- Python 3.11+ with eFISHent and all Python deps
- bowtie2, jellyfish, blastn, dustmasker (from BLAST+), glpsol
- Fold binary (Linux)
- Entrypoint script that: downloads indices from S3 → runs eFISHent → uploads results to S3

Image size estimate: ~1–2 GB (stored in ECR).

---

## Pricing Breakdown (us-east-1)

### Fixed Monthly Costs (always-on)

| Component | Spec | Cost/month |
|-----------|------|-----------|
| Frontend EC2 | t4g.small (2 vCPU, 2 GB) | **$0** (free tier through Dec 2026), then ~$12/mo |
| EBS volume | 20 GB gp3 (frontend OS + app) | $1.60 |
| S3 storage | ~10 GB (genome indices + results) | $0.23 |
| ECR | ~2 GB Docker image | $0.20 |
| Public IPv4 | 1 address | $3.60 |
| **Total fixed** | | **~$6/mo** (free tier) or **~$18/mo** (after) |

### Per-Job Costs (Spot compute)

| Job Type | Instance | vCPU | RAM | Spot $/hr | Typical Duration | Cost/job |
|----------|----------|------|-----|-----------|-----------------|----------|
| Small gene (<5kb) | c6i.2xlarge | 8 | 16 GB | ~$0.10 | 10–15 min | **$0.02–0.03** |
| Medium gene (5–50kb) | c6i.2xlarge | 8 | 16 GB | ~$0.10 | 20–45 min | **$0.03–0.08** |
| Large gene (>50kb) | c6i.4xlarge | 16 | 32 GB | ~$0.20 | 30–90 min | **$0.10–0.30** |

Spot pricing: c6i instances typically run at ~70–85% discount off on-demand ($0.34/hr → ~$0.07–0.10/hr spot for 2xlarge).

S3 transfer between Batch and S3 in the same region is **free**.

### Monthly Cost Scenarios

| Usage | Jobs/month | Avg job cost | Compute | Fixed | **Total** |
|-------|-----------|-------------|---------|-------|-----------|
| Light (personal/lab) | 10 | $0.05 | $0.50 | $6 | **~$7/mo** |
| Moderate (small group) | 50 | $0.05 | $2.50 | $6 | **~$9/mo** |
| Active (course/workshop) | 200 | $0.08 | $16.00 | $6 | **~$22/mo** |
| Heavy (public service) | 1000 | $0.08 | $80.00 | $6 | **~$86/mo** |

> **Bottom line: Running a public eFISHent service for a typical lab costs under $10/month.**

---

## Implementation Steps

### Phase 1: Docker Image

```dockerfile
FROM python:3.11-slim

# Install system dependencies
RUN apt-get update && apt-get install -y \
    bowtie2 jellyfish ncbi-blast+ glpk-utils wget \
    && rm -rf /var/lib/apt/lists/*

# Install eFISHent
COPY . /app
RUN pip install /app

# Copy Fold binary
COPY eFISHent/Fold_linux /usr/local/bin/Fold
RUN chmod +x /usr/local/bin/Fold

# Entrypoint handles S3 sync + run + upload
COPY hosting/entrypoint.sh /entrypoint.sh
ENTRYPOINT ["/entrypoint.sh"]
```

### Phase 2: Pre-Build Genome Indices

```bash
# On a c6i.4xlarge (one-time, ~30 min)
eFISHent --reference-genome Homo_sapiens.GRCh38.dna.primary_assembly.fa --build-indices True
eFISHent --reference-genome Mus_musculus.GRCm39.dna.primary_assembly.fa --build-indices True

# Upload to S3
aws s3 sync ./indices/ s3://efishent-genomes/indices/
```

### Phase 3: AWS Batch Setup

1. **ECR repository** — push Docker image
2. **Job definition** — container image, vCPU/memory requirements, environment variables (S3 paths, job ID)
3. **Compute environment** — Spot, c6i family, min 0 / max 4 instances, `SPOT_PRICE_CAPACITY_OPTIMIZED` allocation strategy
4. **Job queue** — single queue, priority 1

### Phase 4: Frontend

Minimal FastAPI app:

- `GET /` — form: gene name or FASTA upload, organism dropdown, preset dropdown, advanced params accordion
- `POST /submit` — validate, create job record, call `boto3.client('batch').submit_job()`
- `GET /status/{job_id}` — poll Batch job status, show progress
- `GET /results/{job_id}` — pre-signed S3 download links for output files
- Optional: email notification via SES when job completes ($0.10 per 1000 emails)

### Phase 5: Cleanup & Limits

- S3 lifecycle rule: delete result files after 7 days
- Job timeout: 2 hours max
- Rate limit: 5 concurrent jobs, 20 jobs/day per IP (prevent abuse)
- Batch max instances: 4 (caps spend at ~$0.80/hr worst case)

---

## What's NOT Needed for MVP

- **No database server** — SQLite on the frontend EBS or DynamoDB (25 GB free tier) for job metadata
- **No Redis/Celery** — AWS Batch IS the queue
- **No Kubernetes** — massive overkill for this workload
- **No user accounts** — job IDs are enough for MVP; add auth later if needed
- **No custom domain/SSL** — use the EC2 public IP or a free Cloudflare tunnel initially; add a domain later

---

## Optional Enhancements (Post-MVP)

| Enhancement | Effort | Value |
|------------|--------|-------|
| Custom domain + HTTPS (via Cloudflare or ACM + ALB) | 1–2 hours | Professional appearance |
| Email notifications (SES) | 2–3 hours | Users don't have to poll |
| Pre-computed probes cache (hash params → cached results) | 3–4 hours | Instant results for popular genes |
| User accounts + job history | 1 day | Returning users see past jobs |
| Additional organisms (Drosophila, zebrafish, C. elegans) | 1 hour each | Broader audience |
| Usage dashboard (CloudWatch) | 1–2 hours | Monitor costs and usage |
| GPU-accelerated alignment (bowtie2 on GPU) | Not worth it | Marginal speedup, high cost |

---

## Frontend Design

Target audience: bench biologists who run FISH experiments. Comfortable with basic bioinformatics but don't want to think about CLI flags. The UI should feel like a modern lab tool, not a developer dashboard.

### Pages

**Landing page** (`/`):
- Clean hero with eFISHent logo, one-liner ("Design smFISH probes in minutes"), prominent "Design Probes" CTA
- Value prop cards: specificity (7-layer filtering), optimized spacing, quality-scored output
- Live example: pre-computed ACTB probe map that users can interact with before submitting anything
- Footer: citation/DOI, GitHub link, documentation

**Job submission** (`/design`) — single page, not a wizard:
- **Target**: Toggle between "Gene name" (autocomplete from Ensembl) and "Upload FASTA" (drag-and-drop). Gene name mode shows suggestions as you type — eliminates the biggest friction point
- **Organism**: Large selectable cards for Human (GRCh38) and Mouse (GRCm39). "Request organism" link opens a GitHub issue template
- **Protocol preset**: Selectable cards for smFISH, MERFISH, DNA-FISH, Custom — each with a short description. Selecting one fills defaults
- **Advanced parameters** (collapsed by default): Accordion sections for probe length/spacing, thermodynamic filters, off-target filters, optimization. Each param has current value, range indicator, plain-language tooltip. Presets grey out overridden params with "reset to preset" link
- **Submit**: Shows estimated runtime ("~15 min for a 3kb gene"). Optional email field for notification
- **Live validation**: Red outline + message for conflicts (min_tm > max_tm, etc.)

**Job status** (`/job/{id}`):
- Progress bar mapped to pipeline stages (Candidates → Basic filtering → Alignment → K-mer → Structure → Optimization → Done)
- Each completed stage shows a stat: "12,450 candidates → 8,230 passed basic → 2,104 passed alignment → ..."
- Estimated time remaining based on gene length + current stage
- Browser notification (Notification API) when done — users can switch tabs

**Results** (`/job/{id}` — same URL, updates when done):
- **Probe map** (star feature): SVG visualization of gene with probes as colored bars. Color = quality score (green → yellow → red). Hover shows Tm, GC, quality, sequence. Zoomable/pannable
- **Summary cards**: Probe count, avg quality, coverage %, avg spacing, Tm range
- **Probe table**: Sortable by quality/Tm/GC/position. Click row to highlight on map. Checkbox to include/exclude individual probes
- **Downloads**: FASTA, CSV, PDF report, and "Copy for IDT order" button (formats sequences for IDT bulk ordering — huge time-saver)
- **Share**: Copy job URL. Results persist 7 days

### Visual Style

- Clean, generous white space, subtle shadows — Linear/Vercel design language
- Accent color from eFISHent logo palette
- System font stack (no custom fonts to load)
- Dark mode toggle (microscopists work in dark rooms)
- No gratuitous animations; smooth transitions between states
- Mobile-friendly (biologists check results on phone between experiments)

### Tech Stack

- **Option A (simple)**: HTMX + Jinja templates — no build step, server-rendered, fast. Best for MVP
- **Option B (richer)**: React + Tailwind — better for the interactive probe map. Adds build complexity
- **Probe map**: D3.js or lightweight Canvas library
- No component library — keep bundle small, design custom

### Key Differentiators

The probe map + IDT order formatting are the two features that would make this genuinely more useful than wrapping the CLI in a form. Pre-computed example results (ACTB, GAPDH) let users see what they'll get before waiting.

---

## Risk Considerations

| Risk | Mitigation |
|------|-----------|
| Spot interruption mid-job | Batch auto-retries (set `attempts: 3`). Jobs are idempotent |
| Abuse / crypto mining | Container only runs eFISHent entrypoint. Rate limits. Max instance cap |
| Large gene takes >2 hours | Timeout + notify user. Suggest splitting or using `--optimization-method greedy` |
| S3 costs from large results | Lifecycle policy deletes after 7 days. Results are small (~10 MB) |
| Genome index download slow at job start | Use EFS instead of S3 for indices if latency matters (adds ~$5/mo for 10 GB) |
