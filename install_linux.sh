#!/bin/bash
# eFISHent Linux Installation Script (No sudo required)
# Installs all dependencies to a custom path or ~/.local/efishent-deps
#
# Usage: ./install_linux.sh [INSTALL_DIR]
#   INSTALL_DIR: Optional custom installation directory (default: ~/.local/efishent-deps)

set -e

# Configuration - use first argument or default
INSTALL_DIR="${1:-${HOME}/.local/efishent-deps}"
INSTALL_DIR=$(cd "$(dirname "${INSTALL_DIR}")" 2>/dev/null && pwd)/$(basename "${INSTALL_DIR}") || INSTALL_DIR=$(mkdir -p "${INSTALL_DIR}" && cd "${INSTALL_DIR}" && pwd)
BIN_DIR="${INSTALL_DIR}/bin"
EDIRECT_DIR="${INSTALL_DIR}/edirect"
mkdir -p "${BIN_DIR}"

echo "Installing eFISHent dependencies to ${INSTALL_DIR}"
echo "=================================================="

# Detect architecture
ARCH=$(uname -m)
if [ "$ARCH" != "x86_64" ] && [ "$ARCH" != "aarch64" ]; then
    echo "Unsupported architecture: $ARCH"
    exit 1
fi

cd "${INSTALL_DIR}"

# 1. Bowtie (pre-compiled)
echo ""
echo "[1/7] Installing Bowtie..."
if [ ! -f "${BIN_DIR}/bowtie" ]; then
    wget -q "https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.3.1/bowtie-1.3.1-linux-x86_64.zip/download" -O bowtie.zip
    unzip -q bowtie.zip
    cp bowtie-1.3.1-linux-x86_64/bowtie* "${BIN_DIR}/"
    rm -rf bowtie.zip bowtie-1.3.1-linux-x86_64
    echo "  Bowtie installed successfully"
else
    echo "  Bowtie already installed"
fi

# 2. Bowtie2 (pre-compiled)
echo ""
echo "[2/7] Installing Bowtie2..."
if [ ! -f "${BIN_DIR}/bowtie2" ]; then
    wget -q "https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.5.4/bowtie2-2.5.4-linux-x86_64.zip/download" -O bowtie2.zip
    unzip -q bowtie2.zip
    cp bowtie2-2.5.4-linux-x86_64/bowtie2* "${BIN_DIR}/"
    rm -rf bowtie2.zip bowtie2-2.5.4-linux-x86_64
    echo "  Bowtie2 installed successfully"
else
    echo "  Bowtie2 already installed"
fi

# 3. Jellyfish (compile from source)
echo ""
echo "[3/7] Installing Jellyfish..."
if [ ! -f "${BIN_DIR}/jellyfish" ]; then
    echo "  Downloading and compiling (this may take a minute)..."
    wget -q "https://github.com/gmarcais/Jellyfish/releases/download/v2.3.1/jellyfish-2.3.1.tar.gz"
    tar -xzf jellyfish-2.3.1.tar.gz
    cd jellyfish-2.3.1
    ./configure --prefix="${INSTALL_DIR}" >/dev/null 2>&1
    make -j$(nproc) >/dev/null 2>&1
    make install >/dev/null 2>&1
    cd "${INSTALL_DIR}"
    rm -rf jellyfish-2.3.1.tar.gz jellyfish-2.3.1
    echo "  Jellyfish installed successfully"
else
    echo "  Jellyfish already installed"
fi

# 4. GLPK (compile from source - needed for pyomo solver)
echo ""
echo "[4/7] Installing GLPK..."
if [ ! -f "${BIN_DIR}/glpsol" ]; then
    echo "  Downloading and compiling (this may take a minute)..."
    wget -q "https://ftp.gnu.org/gnu/glpk/glpk-5.0.tar.gz"
    tar -xzf glpk-5.0.tar.gz
    cd glpk-5.0
    ./configure --prefix="${INSTALL_DIR}" >/dev/null 2>&1
    make -j$(nproc) >/dev/null 2>&1
    make install >/dev/null 2>&1
    cd "${INSTALL_DIR}"
    rm -rf glpk-5.0.tar.gz glpk-5.0
    echo "  GLPK installed successfully"
else
    echo "  GLPK already installed"
fi

# 5. BLAST+ (pre-compiled)
echo ""
echo "[5/7] Installing BLAST+..."
if [ ! -f "${BIN_DIR}/blastn" ]; then
    wget -q "https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.16.0+-x64-linux.tar.gz" -O blast.tar.gz
    tar -xzf blast.tar.gz
    cp ncbi-blast-*/bin/blastn ncbi-blast-*/bin/makeblastdb "${BIN_DIR}/"
    rm -rf blast.tar.gz ncbi-blast-*
    echo "  BLAST+ installed successfully"
else
    echo "  BLAST+ already installed"
fi

# 6. Entrez Direct (official installer)
echo ""
echo "[6/7] Installing Entrez Direct..."
if [ ! -f "${EDIRECT_DIR}/esearch" ]; then
    echo "  Downloading and setting up..."
    mkdir -p "${EDIRECT_DIR}"
    wget -q "https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/edirect.tar.gz"
    tar -xzf edirect.tar.gz -C "${EDIRECT_DIR}" --strip-components=1
    rm edirect.tar.gz
    # Run setup silently (it downloads additional components)
    "${EDIRECT_DIR}/setup.sh" >/dev/null 2>&1 || true
    echo "  Entrez Direct installed successfully"
else
    echo "  Entrez Direct already installed"
fi

# 7. RNAstructure Fold (Linux only - macOS has bundled binary)
echo ""
echo "[7/7] Installing RNAstructure Fold..."
if [ ! -f "${BIN_DIR}/Fold" ]; then
    echo "  RNAstructure requires manual download due to registration requirement."
    echo "  Please:"
    echo "    1. Visit: https://rna.urmc.rochester.edu/RNAstructure.html"
    echo "    2. Download 'Text interfaces for 64-bit Linux'"
    echo "    3. Extract and copy 'exe/Fold' to ${BIN_DIR}/"
    echo ""
    echo "  Alternatively, you can try the direct link (may require registration):"
    echo "    wget 'https://rna.urmc.rochester.edu/RNAstructureDownload/RNAstructure_Linux64.tar.gz'"
    echo "    tar -xzf RNAstructure_Linux64.tar.gz"
    echo "    cp RNAstructure/exe/Fold ${BIN_DIR}/"
else
    echo "  RNAstructure Fold already installed"
fi

# Create activation script with actual paths
cat > "${INSTALL_DIR}/activate.sh" << EOF
# Source this file to add eFISHent dependencies to PATH
# Usage: source ${INSTALL_DIR}/activate.sh

export PATH="${BIN_DIR}:${EDIRECT_DIR}:\${PATH}"
export LD_LIBRARY_PATH="${INSTALL_DIR}/lib:\${LD_LIBRARY_PATH}"
EOF

# Activate and verify installation
# shellcheck source=/dev/null
source "${INSTALL_DIR}/activate.sh"

echo ""
echo "=================================================="
echo "Installation complete! Installed versions:"
echo ""
bowtie_ver="$(bowtie --version 2>&1 | head -1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+')"
bowtie2_ver="$(bowtie2 --version 2>&1 | head -1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+')"
jellyfish_ver="$(jellyfish --version 2>&1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+')"
blastn_ver="$(blastn -version 2>&1 | head -1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+')"
glpsol_ver="$(glpsol --version 2>&1 | head -1 | grep -oE '[0-9]+\.[0-9]+')"
edirect_ver="$(esearch -version 2>&1 | head -1)"
echo "  bowtie:      $bowtie_ver"
echo "  bowtie2:     $bowtie2_ver"
echo "  jellyfish:   $jellyfish_ver"
echo "  blastn:      $blastn_ver"
echo "  glpsol:      $glpsol_ver"
echo "  edirect:     $edirect_ver"
if [ -f "${BIN_DIR}/Fold" ]; then
    echo "  Fold:        installed"
else
    echo "  Fold:        NOT INSTALLED (see instructions above)"
fi
echo ""
echo "To use eFISHent, add to your ~/.bashrc:"
echo ""
echo "  source ${INSTALL_DIR}/activate.sh"
echo ""
echo "Then install the Python package:"
echo ""
echo "  pip install efishent"
echo ""
echo "Or using uv (faster alternative):"
echo ""
echo "  uv pip install efishent"
echo ""
echo "  # Or create a managed virtual environment:"
echo "  uv venv && source .venv/bin/activate && uv pip install efishent"
echo ""
