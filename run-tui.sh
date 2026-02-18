#!/bin/bash
set -euo pipefail

# Run Hedgehog TUI
cd "$(dirname "$0")/tui"

# Ensure TypeScript compiler is available from local devDependencies.
# This typically fails when dependencies were installed in production mode.
if [ ! -x "node_modules/.bin/tsc" ]; then
  echo "TypeScript compiler not found. Installing TUI dependencies..."
  npm_config_production=false npm ci
fi

npm run build
npm start
