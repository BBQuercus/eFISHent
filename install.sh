#!/bin/bash
# eFISHent — One-Command Installer
# Installs eFISHent and all dependencies. No sudo, Docker, or conda required.
#
# Usage:
#   curl -LsSf https://raw.githubusercontent.com/BBQuercus/eFISHent/main/install.sh | sh
#   curl -LsSf https://raw.githubusercontent.com/BBQuercus/eFISHent/main/install.sh | sh -s -- --prefix ~/.local/efishent
#
# Options:
#   --prefix DIR    Install to DIR (default: ~/.local/efishent)
#   --no-modify-rc  Don't modify shell rc files
#   --deps-only     Only install external dependencies (skip Python package)
#   --uninstall     Remove eFISHent and all dependencies

set -e

# ── Configuration ─────────────────────────────────────────────────────────────

EFISHENT_VERSION="latest"
BOWTIE_VERSION="1.3.1"
JELLYFISH_VERSION="2.3.1"
GLPK_VERSION="5.0"

# Colors (disabled if not a terminal)
if [ -t 1 ]; then
    BOLD=$(printf '\033[1m')
    DIM=$(printf '\033[2m')
    GREEN=$(printf '\033[32m')
    YELLOW=$(printf '\033[33m')
    RED=$(printf '\033[31m')
    CYAN=$(printf '\033[36m')
    RESET=$(printf '\033[0m')
else
    BOLD="" DIM="" GREEN="" YELLOW="" RED="" CYAN="" RESET=""
fi

# ── Helpers ───────────────────────────────────────────────────────────────────

info()  { printf "%s %s\n" "${GREEN}✓${RESET}" "$1"; }
warn()  { printf "%s %s\n" "${YELLOW}!${RESET}" "$1"; }
err()   { printf "%s %s\n" "${RED}✗${RESET}" "$1" >&2; }
step()  { printf "\n%s\n" "${BOLD}$1${RESET}"; }

need_cmd() {
    if ! command -v "$1" >/dev/null 2>&1; then
        err "Required command not found: $1"
        err "Please install it and retry."
        exit 1
    fi
}

# Download helper — uses curl or wget
download() {
    local url="$1" dest="$2"
    if command -v curl >/dev/null 2>&1; then
        curl -fsSL "$url" -o "$dest"
    elif command -v wget >/dev/null 2>&1; then
        wget -q "$url" -O "$dest"
    else
        err "Neither curl nor wget found. Install one and retry."
        exit 1
    fi
}

# ── Parse arguments ───────────────────────────────────────────────────────────

INSTALL_PREFIX="${HOME}/.local/efishent"
MODIFY_RC=true
DEPS_ONLY=false
UNINSTALL=false

while [ $# -gt 0 ]; do
    case "$1" in
        --prefix)    INSTALL_PREFIX="$2"; shift 2 ;;
        --no-modify-rc) MODIFY_RC=false; shift ;;
        --deps-only) DEPS_ONLY=true; shift ;;
        --uninstall) UNINSTALL=true; shift ;;
        *)           err "Unknown option: $1"; exit 1 ;;
    esac
done

DEPS_DIR="${INSTALL_PREFIX}/deps"
BIN_DIR="${DEPS_DIR}/bin"
VENV_DIR="${INSTALL_PREFIX}/venv"
WRAPPER_DIR="${HOME}/.local/bin"
EDIRECT_DIR="${DEPS_DIR}/edirect"

# ── Uninstall ─────────────────────────────────────────────────────────────────

if [ "$UNINSTALL" = true ]; then
    step "Uninstalling eFISHent..."
    rm -rf "${INSTALL_PREFIX}"
    rm -f "${WRAPPER_DIR}/efishent"
    info "Removed ${INSTALL_PREFIX}"
    info "Removed ${WRAPPER_DIR}/efishent"
    warn "You may want to remove the PATH entry from your shell rc file."
    exit 0
fi

# ── Detect platform ──────────────────────────────────────────────────────────

OS=$(uname -s | tr '[:upper:]' '[:lower:]')
ARCH=$(uname -m)

case "$OS" in
    linux)  PLATFORM="linux" ;;
    darwin) PLATFORM="macos" ;;
    *)      err "Unsupported OS: $OS"; exit 1 ;;
esac

case "$ARCH" in
    x86_64)         ARCH_NAME="x86_64" ;;
    aarch64|arm64)  ARCH_NAME="aarch64" ;;
    *)              err "Unsupported architecture: $ARCH"; exit 1 ;;
esac

# ── Banner ────────────────────────────────────────────────────────────────────

printf "\n%s\n" "${BOLD}eFISHent Installer${RESET}"
printf "%s\n" "${DIM}Platform: ${PLATFORM}/${ARCH_NAME}${RESET}"
printf "%s\n" "${DIM}Install:  ${INSTALL_PREFIX}${RESET}"

# ── Prerequisites ─────────────────────────────────────────────────────────────

# Need at least curl or wget
if ! command -v curl >/dev/null 2>&1 && ! command -v wget >/dev/null 2>&1; then
    err "Neither curl nor wget found. Install one and retry."
    exit 1
fi

# Need tar for extracting archives
need_cmd tar

# Need unzip for bowtie
need_cmd unzip

# Check for a C compiler (still needed for jellyfish/glpk if no binary available)
HAS_COMPILER=false
if command -v gcc >/dev/null 2>&1 || command -v clang >/dev/null 2>&1; then
    HAS_COMPILER=true
fi
if command -v make >/dev/null 2>&1; then
    HAS_MAKE=true
else
    HAS_MAKE=false
fi

mkdir -p "${BIN_DIR}" "${DEPS_DIR}" "${WRAPPER_DIR}"

# ── 1. Install Bowtie ────────────────────────────────────────────────────────

step "[1/5] Bowtie (sequence aligner)"

if [ -f "${BIN_DIR}/bowtie" ]; then
    info "Already installed"
else
    case "${PLATFORM}-${ARCH_NAME}" in
        linux-x86_64)
            BOWTIE_URL="https://sourceforge.net/projects/bowtie-bio/files/bowtie/${BOWTIE_VERSION}/bowtie-${BOWTIE_VERSION}-linux-x86_64.zip/download"
            BOWTIE_DIR="bowtie-${BOWTIE_VERSION}-linux-x86_64"
            ;;
        macos-x86_64)
            BOWTIE_URL="https://sourceforge.net/projects/bowtie-bio/files/bowtie/${BOWTIE_VERSION}/bowtie-${BOWTIE_VERSION}-macos-x86_64.zip/download"
            BOWTIE_DIR="bowtie-${BOWTIE_VERSION}-macos-x86_64"
            ;;
        macos-aarch64)
            # Bowtie 1 doesn't have native ARM builds — x86_64 works via Rosetta
            BOWTIE_URL="https://sourceforge.net/projects/bowtie-bio/files/bowtie/${BOWTIE_VERSION}/bowtie-${BOWTIE_VERSION}-macos-x86_64.zip/download"
            BOWTIE_DIR="bowtie-${BOWTIE_VERSION}-macos-x86_64"
            warn "Using x86_64 build (runs via Rosetta 2 on Apple Silicon)"
            ;;
        linux-aarch64)
            # No pre-built aarch64 bowtie — check if it's already on PATH
            if command -v bowtie >/dev/null 2>&1; then
                info "Using system bowtie: $(command -v bowtie)"
                ln -sf "$(command -v bowtie)" "${BIN_DIR}/bowtie"
                BOWTIE_URL=""
            else
                err "No pre-built Bowtie for Linux aarch64."
                err "Install bowtie via your package manager and retry."
                exit 1
            fi
            ;;
    esac

    if [ -n "${BOWTIE_URL:-}" ]; then
        download "$BOWTIE_URL" "${DEPS_DIR}/bowtie.zip"
        cd "${DEPS_DIR}"
        unzip -qo bowtie.zip
        cp "${BOWTIE_DIR}"/bowtie* "${BIN_DIR}/"
        rm -rf bowtie.zip "${BOWTIE_DIR}"
        # Remove macOS quarantine
        if [ "$PLATFORM" = "macos" ]; then
            xattr -dr com.apple.quarantine "${BIN_DIR}/bowtie"* 2>/dev/null || true
        fi
        cd - >/dev/null
        info "Installed bowtie ${BOWTIE_VERSION}"
    fi
fi

# ── 2. Install Jellyfish ─────────────────────────────────────────────────────

step "[2/5] Jellyfish (k-mer counter)"

if [ -f "${BIN_DIR}/jellyfish" ]; then
    info "Already installed"
else
    # Try pre-compiled binary from conda-forge first (static builds)
    JELLYFISH_INSTALLED=false

    # Fallback: compile from source
    if [ "$JELLYFISH_INSTALLED" = false ]; then
        if [ "$HAS_COMPILER" = true ] && [ "$HAS_MAKE" = true ]; then
            printf "  Compiling from source (this may take a minute)...\n"
            download "https://github.com/gmarcais/Jellyfish/releases/download/v${JELLYFISH_VERSION}/jellyfish-${JELLYFISH_VERSION}.tar.gz" "${DEPS_DIR}/jellyfish.tar.gz"
            cd "${DEPS_DIR}"
            tar -xzf jellyfish.tar.gz
            cd "jellyfish-${JELLYFISH_VERSION}"

            NPROC=1
            if command -v nproc >/dev/null 2>&1; then
                NPROC=$(nproc)
            elif command -v sysctl >/dev/null 2>&1; then
                NPROC=$(sysctl -n hw.ncpu 2>/dev/null || echo 1)
            fi

            if ./configure --prefix="${DEPS_DIR}" > "${DEPS_DIR}/jellyfish_build.log" 2>&1 && \
               make -j"${NPROC}" >> "${DEPS_DIR}/jellyfish_build.log" 2>&1 && \
               make install >> "${DEPS_DIR}/jellyfish_build.log" 2>&1; then
                info "Compiled and installed jellyfish ${JELLYFISH_VERSION}"
                JELLYFISH_INSTALLED=true
            else
                err "Jellyfish compilation failed. See ${DEPS_DIR}/jellyfish_build.log"
                exit 1
            fi
            cd "${DEPS_DIR}"
            rm -rf jellyfish.tar.gz "jellyfish-${JELLYFISH_VERSION}"
        else
            err "No C compiler or make found — cannot build jellyfish."
            err "Install gcc and make, or install jellyfish manually."
            exit 1
        fi
    fi
fi

# ── 3. Install GLPK ──────────────────────────────────────────────────────────

step "[3/5] GLPK (linear programming solver)"

if [ -f "${BIN_DIR}/glpsol" ]; then
    info "Already installed"
elif command -v glpsol >/dev/null 2>&1; then
    info "Using system glpsol: $(command -v glpsol)"
    ln -sf "$(command -v glpsol)" "${BIN_DIR}/glpsol"
else
    if [ "$HAS_COMPILER" = true ] && [ "$HAS_MAKE" = true ]; then
        printf "  Compiling from source (this may take a minute)...\n"
        download "https://ftp.gnu.org/gnu/glpk/glpk-${GLPK_VERSION}.tar.gz" "${DEPS_DIR}/glpk.tar.gz"
        cd "${DEPS_DIR}"
        tar -xzf glpk.tar.gz
        cd "glpk-${GLPK_VERSION}"

        NPROC=1
        if command -v nproc >/dev/null 2>&1; then
            NPROC=$(nproc)
        elif command -v sysctl >/dev/null 2>&1; then
            NPROC=$(sysctl -n hw.ncpu 2>/dev/null || echo 1)
        fi

        if ./configure --prefix="${DEPS_DIR}" > "${DEPS_DIR}/glpk_build.log" 2>&1 && \
           make -j"${NPROC}" >> "${DEPS_DIR}/glpk_build.log" 2>&1 && \
           make install >> "${DEPS_DIR}/glpk_build.log" 2>&1; then
            info "Compiled and installed glpsol ${GLPK_VERSION}"
        else
            err "GLPK compilation failed. See ${DEPS_DIR}/glpk_build.log"
            exit 1
        fi
        cd "${DEPS_DIR}"
        rm -rf glpk.tar.gz "glpk-${GLPK_VERSION}"
    else
        warn "No compiler found — skipping GLPK."
        warn "glpsol is only needed for --optimization-method optimal (greedy works without it)."
    fi
fi

# ── 4. Install Entrez Direct ─────────────────────────────────────────────────

step "[4/5] Entrez Direct (NCBI tools — optional)"

if [ -f "${EDIRECT_DIR}/esearch" ]; then
    info "Already installed"
elif command -v esearch >/dev/null 2>&1; then
    info "Using system esearch: $(command -v esearch)"
else
    # Backup existing ~/edirect if present
    EDIRECT_BACKUP=""
    if [ -d "${HOME}/edirect" ]; then
        EDIRECT_BACKUP="${HOME}/edirect.bak.$$"
        mv "${HOME}/edirect" "$EDIRECT_BACKUP"
    fi

    INSTALL_CMD=""
    if command -v curl >/dev/null 2>&1; then
        INSTALL_CMD='sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"'
    elif command -v wget >/dev/null 2>&1; then
        INSTALL_CMD='sh -c "$(wget -q https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh -O -)"'
    fi

    if eval "$INSTALL_CMD" > "${DEPS_DIR}/edirect_install.log" 2>&1; then
        if [ -d "${HOME}/edirect" ]; then
            mv "${HOME}/edirect" "${EDIRECT_DIR}"
            info "Installed Entrez Direct"
        fi
    else
        warn "Entrez Direct installation failed (optional — only needed for NCBI gene download)."
        warn "You can provide sequence files directly with --sequence-file instead."
    fi

    # Restore backup
    if [ -n "$EDIRECT_BACKUP" ] && [ -d "$EDIRECT_BACKUP" ]; then
        mv "$EDIRECT_BACKUP" "${HOME}/edirect"
    fi
fi

# ── 5. Install eFISHent Python package ────────────────────────────────────────

if [ "$DEPS_ONLY" = false ]; then
    step "[5/5] eFISHent (Python package)"

    # Install uv if not present
    if ! command -v uv >/dev/null 2>&1; then
        printf "  Installing uv (Python package manager)...\n"
        if command -v curl >/dev/null 2>&1; then
            curl -LsSf https://astral.sh/uv/install.sh | sh 2>/dev/null
        else
            wget -qO- https://astral.sh/uv/install.sh | sh 2>/dev/null
        fi
        # Source cargo env to get uv on PATH
        [ -f "${HOME}/.local/bin/uv" ] && export PATH="${HOME}/.local/bin:${PATH}"
        [ -f "${HOME}/.cargo/env" ] && . "${HOME}/.cargo/env"
        if ! command -v uv >/dev/null 2>&1; then
            err "uv installation failed. Install manually: https://docs.astral.sh/uv/"
            exit 1
        fi
        info "Installed uv"
    fi

    # Create venv and install efishent
    if [ ! -d "${VENV_DIR}" ]; then
        uv venv "${VENV_DIR}" --python ">=3.10" > /dev/null 2>&1
        info "Created Python virtual environment"
    fi

    printf "  Installing eFISHent...\n"
    uv pip install --python "${VENV_DIR}/bin/python" efishent 2>/dev/null
    info "Installed eFISHent"
else
    step "[5/5] Skipping Python package (--deps-only)"
fi

# ── Create wrapper script ────────────────────────────────────────────────────

step "Creating efishent wrapper..."

# Determine library path variable
if [ "$PLATFORM" = "macos" ]; then
    LIB_PATH_VAR="DYLD_LIBRARY_PATH"
else
    LIB_PATH_VAR="LD_LIBRARY_PATH"
fi

cat > "${WRAPPER_DIR}/efishent" << WRAPPER
#!/bin/sh
# eFISHent wrapper — automatically sets up PATH for dependencies
# Generated by install.sh — do not edit manually

# Add dependency binaries to PATH
export PATH="${BIN_DIR}:${EDIRECT_DIR}:\${PATH}"
export ${LIB_PATH_VAR}="${DEPS_DIR}/lib:\${${LIB_PATH_VAR}}"

# Run eFISHent from its venv
exec "${VENV_DIR}/bin/efishent" "\$@"
WRAPPER
chmod +x "${WRAPPER_DIR}/efishent"
info "Created wrapper at ${WRAPPER_DIR}/efishent"

# Also create an activate script for users who want to use the deps in their shell
cat > "${INSTALL_PREFIX}/activate.sh" << ACTIVATE
# Source this to add eFISHent dependencies to your current shell
# Usage: source ${INSTALL_PREFIX}/activate.sh

export PATH="${BIN_DIR}:${EDIRECT_DIR}:\${PATH}"
export ${LIB_PATH_VAR}="${DEPS_DIR}/lib:\${${LIB_PATH_VAR}}"

# Also make the venv's efishent available
export PATH="${VENV_DIR}/bin:\${PATH}"
ACTIVATE
info "Created activate script at ${INSTALL_PREFIX}/activate.sh"

# ── Add to PATH ──────────────────────────────────────────────────────────────

RC_UPDATED=false

if [ "$MODIFY_RC" = true ]; then
    # Determine shell rc file
    SHELL_NAME=$(basename "${SHELL:-/bin/sh}")
    case "$SHELL_NAME" in
        zsh)  RC_FILE="${HOME}/.zshrc" ;;
        bash) RC_FILE="${HOME}/.bashrc" ;;
        fish) RC_FILE="${HOME}/.config/fish/config.fish" ;;
        *)    RC_FILE="${HOME}/.profile" ;;
    esac

    MARKER='# Added by eFISHent installer'

    if [ -f "$RC_FILE" ] && grep -q "$MARKER" "$RC_FILE" 2>/dev/null; then
        info "PATH already configured in ${RC_FILE}"
    else
        printf "\n"
        printf "Add %s to PATH in %s? [Y/n] " "$WRAPPER_DIR" "$RC_FILE"

        # Handle non-interactive (piped) input — default to yes
        if [ -t 0 ]; then
            read -r response
        else
            response="y"
            printf "y (non-interactive)\n"
        fi

        case "$response" in
            [nN]*)
                warn "Skipped. Add manually: export PATH=\"${WRAPPER_DIR}:\$PATH\""
                ;;
            *)
                if [ "$SHELL_NAME" = "fish" ]; then
                    mkdir -p "$(dirname "$RC_FILE")"
                    printf '\n%s\nfish_add_path "%s"\n' "$MARKER" "${WRAPPER_DIR}" >> "$RC_FILE"
                else
                    printf '\n%s\nexport PATH="%s:$PATH"\n' "$MARKER" "${WRAPPER_DIR}" >> "$RC_FILE"
                fi
                info "Added to ${RC_FILE}"
                RC_UPDATED=true
                ;;
        esac
    fi
fi

# ── Verify installation ──────────────────────────────────────────────────────

step "Verifying installation..."

# Temporarily set PATH for verification
export PATH="${BIN_DIR}:${EDIRECT_DIR}:${WRAPPER_DIR}:${PATH}"
export ${LIB_PATH_VAR}="${DEPS_DIR}/lib:$(eval echo "\${${LIB_PATH_VAR}}")"

CHECKS_PASSED=true

for tool in bowtie jellyfish; do
    if command -v "$tool" >/dev/null 2>&1; then
        info "$tool: $(command -v "$tool")"
    else
        err "$tool: not found"
        CHECKS_PASSED=false
    fi
done

# Optional tools
for tool in glpsol esearch; do
    if command -v "$tool" >/dev/null 2>&1; then
        info "$tool: $(command -v "$tool")"
    else
        warn "$tool: not found (optional)"
    fi
done

if [ "$DEPS_ONLY" = false ]; then
    if "${VENV_DIR}/bin/efishent" --version >/dev/null 2>&1; then
        EFISHENT_VER=$("${VENV_DIR}/bin/efishent" --version 2>&1)
        info "efishent: ${EFISHENT_VER}"
    else
        err "efishent: failed to run"
        CHECKS_PASSED=false
    fi
fi

# ── Done ──────────────────────────────────────────────────────────────────────

printf "\n%s\n\n" "${BOLD}${GREEN}Installation complete!${RESET}"

if [ "$RC_UPDATED" = true ]; then
    printf "  Your PATH has been updated in ${CYAN}${RC_FILE}${RESET}.\n"
    printf "  To use efishent in ${BOLD}this${RESET} terminal, run:\n\n"
    printf "    ${CYAN}source ${RC_FILE}${RESET}\n\n"
    printf "  New terminal windows will work automatically.\n\n"
fi

printf "  Get started:\n"
printf "    ${CYAN}efishent --check${RESET}       Verify all dependencies\n"
printf "    ${CYAN}efishent --help${RESET}        Show usage\n"
printf "    ${CYAN}efishent --preset list${RESET}  Show parameter presets\n"
printf "\n"
printf "  To uninstall:\n"
printf "    ${CYAN}curl -LsSf https://raw.githubusercontent.com/BBQuercus/eFISHent/main/install.sh | sh -s -- --uninstall${RESET}\n"
printf "\n"
