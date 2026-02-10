#!/usr/bin/env bash

# SJARACNe_run.sh
# Usage:
#   ./SJARACNe_run.sh [--project_dir DIR] [--project_name NAME]
# Defaults:
#   --project_dir="../"
#   --project_name="NetBID2_Project"

# -------------------------
# Default parameters
# -------------------------
PROJECT_DIR="./"
PROJECT_NAME="NetBID2_Project"

# -------------------------
# Parse command-line options
# -------------------------
while [[ $# -gt 0 ]]; do
    case $1 in
        --project_dir)
            PROJECT_DIR="$2"
            shift 2
            ;;
        --project_name)
            PROJECT_NAME="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            echo "Usage: $0 [--project_dir DIR] [--project_name NAME]"
            exit 1
            ;;
    esac
done

# Ensure PROJECT_DIR ends with a slash
[[ "${PROJECT_DIR}" != */ ]] && PROJECT_DIR="${PROJECT_DIR}/"

# -------------------------
# Construct paths
# -------------------------
BASE_DIR="${PROJECT_DIR}${PROJECT_NAME}/SJAR/${PROJECT_NAME}/"

INPUT_EXP="${BASE_DIR}input.exp"
SIG_FILE="${BASE_DIR}sig.txt"
TF_FILE="${BASE_DIR}tf.txt"
LNC_FILE="${BASE_DIR}lnc_list.txt"

# Intermediate (top 5 only)
SIG_TOP5="${BASE_DIR}sig_top5.txt"
TF_TOP5="${BASE_DIR}tf_top5.txt"
LNC_TOP5="${BASE_DIR}lnc_top5.txt"

# Cleaned files
SIG_CLEAN="${BASE_DIR}sig_cleaned.txt"
TF_CLEAN="${BASE_DIR}tf_cleaned.txt"
LNC_CLEAN="${BASE_DIR}lnc_cleaned.txt"

echo "Using input.exp: $INPUT_EXP"

# Add SJARACNe conda env to PATH
export PATH="/opt/conda/envs/SJARACNe/bin:$PATH"

# Check if sjaracne exists
if ! command -v sjaracne &> /dev/null; then
    echo "Error: sjaracne command not found. Please install or add it to your PATH."
    exit 1
fi

# -------------------------
# Keep top 5 lines
# -------------------------
echo "Keeping top 5 lines of sig.txt, tf.txt, and lnc_list.txt..."
head -n 5 "$SIG_FILE" > "$SIG_TOP5"
head -n 5 "$TF_FILE"  > "$TF_TOP5"
head -n 5 "$LNC_FILE" > "$LNC_TOP5"

# -------------------------
# Clean trailing spaces
# -------------------------
echo "Cleaning trailing spaces..."
sed 's/[[:space:]]*$//' "$SIG_TOP5" > "$SIG_CLEAN"
sed 's/[[:space:]]*$//' "$TF_TOP5"  > "$TF_CLEAN"
sed 's/[[:space:]]*$//' "$LNC_TOP5" > "$LNC_CLEAN"

# -------------------------
# Function to run sjaracne
# -------------------------
run_sjaracne() {
    local cleaned_file="$1"
    local output_dir="$2"

    echo "Running sjaracne for ${cleaned_file} -> ${BASE_DIR}${output_dir}/"

    sjaracne local \
        -e "$INPUT_EXP" \
        -g "$cleaned_file" \
        -o "${BASE_DIR}${output_dir}" \
        -tmp ~/tmp/tmp \
        -n 1

    wait
    rm -rf ~/tmp/tmp

    echo "Finished ${output_dir}."
}

# -------------------------
# Sequential execution
# -------------------------
run_sjaracne "$SIG_CLEAN" "sig"
run_sjaracne "$TF_CLEAN"  "tf"
run_sjaracne "$LNC_CLEAN" "lnc"

echo "All SJARACNe runs completed successfully!"
