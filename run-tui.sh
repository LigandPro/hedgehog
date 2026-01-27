#!/bin/bash
# Run Hedgehog TUI
cd "$(dirname "$0")/tui"
npm run build && npm start
