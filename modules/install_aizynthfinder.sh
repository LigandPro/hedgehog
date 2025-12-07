#!/bin/bash

set -e  

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' 

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
RETROSYNTHESIS_DIR="$SCRIPT_DIR/retrosynthesis"
AIZYNTHFINDER_DIR="$RETROSYNTHESIS_DIR/aizynthfinder"
PUBLIC_DIR="$AIZYNTHFINDER_DIR/public"
DATA_DIR="$AIZYNTHFINDER_DIR/aizynthfinder/data"
LOGGING_SOURCE="$PROJECT_ROOT/src/hedgehog/stages/synthesis/logging.yml"

# Step 1: Check if retrosynthesis repository exists
if [ ! -d "$RETROSYNTHESIS_DIR" ]; then
    echo -e "${YELLOW}  Retrosynthesis repository not found. Cloning...${NC}"
    
    if ! command -v git &> /dev/null; then
        echo -e "${RED} Error: git is not installed. Please install git first.${NC}"
        exit 1
    fi
    
    cd "$SCRIPT_DIR"
    git clone git@github.com:LigandPro/retrosynthesis.git
    echo -e "${GREEN} Repository cloned successfully${NC}"
else
    echo -e "${GREEN} Retrosynthesis repository found${NC}"
fi

# Step 2: Check if aizynthfinder directory exists
if [ ! -d "$AIZYNTHFINDER_DIR" ]; then
    echo -e "${RED} Error: aizynthfinder directory not found in $RETROSYNTHESIS_DIR${NC}"
    exit 1
fi

cd "$AIZYNTHFINDER_DIR"

# Step 3: Check if uv is installed
if ! command -v uv &> /dev/null; then
    echo -e "${YELLOW}  uv is not installed. Installing uv...${NC}"
    curl -LsSf https://astral.sh/uv/install.sh | sh
    export PATH="$HOME/.cargo/bin:$PATH"
    echo -e "${GREEN}✓ uv installed successfully${NC}"
else
    echo -e "${GREEN}✓ uv is installed${NC}"
fi

# Step 4: Install dependencies with uv
echo ""
echo -e "${GREEN}Installing dependencies with uv sync...${NC}"
if uv sync; then
    echo -e "${GREEN} Dependencies installed successfully${NC}"
else
    echo -e "${RED} Error: Failed to install dependencies${NC}"
    exit 1
fi

# Step 5: Download public data 
echo ""
echo -e "${GREEN}Downloading public data ...${NC}"

if [ -d "$PUBLIC_DIR" ] && [ "$(ls -A $PUBLIC_DIR 2>/dev/null)" ]; then
    echo -e "${YELLOW}  Public data directory already exists. Skipping download.${NC}"
    echo -e "${GREEN} Using existing public data${NC}"
else
    echo -e "${GREEN}   Creating public directory and downloading data...${NC}"
    mkdir -p "$PUBLIC_DIR"
    uv run python -m aizynthfinder.tools.download_public_data "$PUBLIC_DIR"
    echo -e "${GREEN} Public data downloaded successfully${NC}"
fi

# Step 6: Copy logging.yml file
echo ""
echo -e "${GREEN}Setting up logging configuration...${NC}"

mkdir -p "$DATA_DIR"

if [ -f "$LOGGING_SOURCE" ]; then
    cp "$LOGGING_SOURCE" "$DATA_DIR/logging.yml"
    echo -e "${GREEN} logging.yml copied successfully${NC}"
else
    if [ ! -f "$DATA_DIR/logging.yml" ]; then
        echo -e "${YELLOW} Source logging.yml not found. Creating default...${NC}"
        cat > "$DATA_DIR/logging.yml" << 'EOF'
version: 1
formatters:
  simple:
    format: '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
handlers:
  console:
    class: logging.StreamHandler
    level: INFO
    formatter: simple
    stream: ext://sys.stdout
  file:
    class: logging.FileHandler
    level: DEBUG
    formatter: simple
    filename: aizynthfinder.log
loggers:
  aizynthfinder:
    level: DEBUG
    handlers: [console, file]
    propagate: no
root:
  level: ERROR
  handlers: [console]
EOF
        echo -e "${GREEN}✓ Default logging.yml created${NC}"
    else
        echo -e "${GREEN}✓ logging.yml already exists${NC}"
    fi
fi

echo ""
echo -e "${GREEN}=========================================="
echo -e "✓ AiZynthFinder installation completed!${NC}"
echo -e "${GREEN}==========================================${NC}"
echo ""

