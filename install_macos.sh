#!/bin/bash
# eFISHent macOS Installation Script (No sudo required)
# Installs all dependencies to a custom path or ~/.local/efishent-deps
# Note: RNAstructure Fold is bundled with eFISHent for macOS
#
# Usage: ./install_macos.sh [INSTALL_DIR]
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
if [ "$ARCH" = "x86_64" ]; then
    ARCH_SUFFIX="x64"
elif [ "$ARCH" = "arm64" ]; then
    ARCH_SUFFIX="aarch64"
else
    echo "Unsupported architecture: $ARCH"
    exit 1
fi

cd "${INSTALL_DIR}"

# 1. Bowtie (pre-compiled)
echo ""
echo "[1/4] Installing Bowtie..."
if [ ! -f "${BIN_DIR}/bowtie" ]; then
    curl -sL "https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.3.1/bowtie-1.3.1-macos-x86_64.zip/download" -o bowtie.zip
    unzip -q bowtie.zip
    cp bowtie-1.3.1-macos-x86_64/bowtie* "${BIN_DIR}/"
    rm -rf bowtie.zip bowtie-1.3.1-macos-x86_64
    # Remove quarantine attribute if present
    xattr -dr com.apple.quarantine "${BIN_DIR}/bowtie"* 2>/dev/null || true
    echo "  Bowtie installed successfully"
else
    echo "  Bowtie already installed"
fi

# 2. Jellyfish (compile from source)
echo ""
echo "[2/4] Installing Jellyfish..."
if [ ! -f "${BIN_DIR}/jellyfish" ]; then
    echo "  Downloading and compiling (this may take a minute)..."
    curl -sL "https://github.com/gmarcais/Jellyfish/releases/download/v2.3.1/jellyfish-2.3.1.tar.gz" -o jellyfish.tar.gz
    tar -xzf jellyfish.tar.gz
    cd jellyfish-2.3.1
    ./configure --prefix="${INSTALL_DIR}" >/dev/null 2>&1
    make -j$(sysctl -n hw.ncpu) >/dev/null 2>&1
    make install >/dev/null 2>&1
    cd "${INSTALL_DIR}"
    rm -rf jellyfish.tar.gz jellyfish-2.3.1
    echo "  Jellyfish installed successfully"
else
    echo "  Jellyfish already installed"
fi

# 3. GLPK (compile from source - needed for pyomo solver)
echo ""
echo "[3/4] Installing GLPK..."
if [ ! -f "${BIN_DIR}/glpsol" ]; then
    echo "  Downloading and compiling (this may take a minute)..."
    curl -sL "https://ftp.gnu.org/gnu/glpk/glpk-5.0.tar.gz" -o glpk.tar.gz
    tar -xzf glpk.tar.gz
    cd glpk-5.0
    ./configure --prefix="${INSTALL_DIR}" >/dev/null 2>&1
    make -j$(sysctl -n hw.ncpu) >/dev/null 2>&1
    make install >/dev/null 2>&1
    cd "${INSTALL_DIR}"
    rm -rf glpk.tar.gz glpk-5.0
    echo "  GLPK installed successfully"
else
    echo "  GLPK already installed"
fi

# 4. Entrez Direct (official installer)
echo ""
echo "[4/4] Installing Entrez Direct..."
if [ ! -f "${EDIRECT_DIR}/esearch" ]; then
    echo "  Downloading and setting up..."
    mkdir -p "${EDIRECT_DIR}"
    curl -sL "https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/edirect.tar.gz" -o edirect.tar.gz
    tar -xzf edirect.tar.gz -C "${EDIRECT_DIR}" --strip-components=1
    rm edirect.tar.gz
    # Run setup silently (it downloads additional components)
    "${EDIRECT_DIR}/setup.sh" >/dev/null 2>&1 || true
    echo "  Entrez Direct installed successfully"
else
    echo "  Entrez Direct already installed"
fi

# Note: RNAstructure Fold is bundled with eFISHent for macOS (Fold_osx)

# Create activation script with actual paths
cat > "${INSTALL_DIR}/activate.sh" << EOF
# Source this file to add eFISHent dependencies to PATH
# Usage: source ${INSTALL_DIR}/activate.sh

export PATH="${BIN_DIR}:${EDIRECT_DIR}:\${PATH}"
export DYLD_LIBRARY_PATH="${INSTALL_DIR}/lib:\${DYLD_LIBRARY_PATH}"
EOF

# Activate and verify installation
source "${INSTALL_DIR}/activate.sh"

echo ""
echo "=================================================="
echo "Installation complete! Installed versions:"
echo ""
echo "  bowtie:      $(bowtie --version 2>&1 | head -1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+')"
echo "  jellyfish:   $(jellyfish --version 2>&1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+')"
echo "  glpsol:      $(glpsol --version 2>&1 | head -1 | grep -oE '[0-9]+\.[0-9]+')"
echo "  edirect:     $(esearch -version 2>&1 | head -1)"
echo "  Fold:        bundled with eFISHent"
echo ""
echo "To use eFISHent, add to your ~/.zshrc (or ~/.bashrc):"
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
