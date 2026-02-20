#!/usr/bin/env bash
# -------------------------
# SJARACNe_run.sh
# -------------------------
# Usage:
#   ./SJARACNe_run.sh [--project_dir DIR] [--project_name NAME]
# Defaults:
#   --project_dir="./"
#   --project_name="NetBID2_Project"
# -------------------------

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
# Input and intermediate files
# -------------------------
INPUT_EXP="${PROJECT_NAME}/input.exp"
SIG_FILE="${PROJECT_NAME}/sig.txt"
TF_FILE="${PROJECT_NAME}/tf.txt"
LNC_FILE="${PROJECT_NAME}/lnc_list.txt"

SIG_TOP5="${PROJECT_NAME}/sig_top5.txt"
TF_TOP5="${PROJECT_NAME}/tf_top5.txt"
LNC_TOP5="${PROJECT_NAME}/lnc_top5.txt"

SIG_CLEAN="${PROJECT_NAME}/sig_cleaned.txt"
TF_CLEAN="${PROJECT_NAME}/tf_cleaned.txt"
LNC_CLEAN="${PROJECT_NAME}/lnc_cleaned.txt"

echo "Using input.exp: $INPUT_EXP"

# -------------------------
# Check input files
# -------------------------
for file in "$SIG_FILE" "$TF_FILE" "$LNC_FILE"; do
    if [ ! -f "$file" ]; then
        echo "Error: Required file not found: $file"
        exit 1
    fi
done

# -------------------------
# Keep top 5 lines
# -------------------------
echo "Keeping top 5 lines of sig.txt, tf.txt, lnc_list.txt"
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
# Add SJARACNe conda env to PATH
# -------------------------
export PATH="/opt/conda/envs/SJARACNe/bin:$PATH"

# Check if sjaracne exists
if ! command -v sjaracne &> /dev/null; then
    echo "Error: sjaracne command not found. Please install or add it to your PATH."
    exit 1
fi

# -------------------------
# Function to run sjaracne
# -------------------------
run_sjaracne() {
    local cleaned_file="$1"
    local output_dir="$2"

    local output_path="${PROJECT_NAME}/${output_dir}"

    mkdir -p "$output_path"

    echo "Running sjaracne for $cleaned_file -> $output_path"

    sjaracne local \
        -e "$INPUT_EXP" \
        -g "$cleaned_file" \
        -o "$output_path" \
        -tmp /tmp/sjaracne/tmp \
        -n 1

    echo "Finished sjaracne run for $output_dir"
}

# -------------------------
# Sequential Execution
# -------------------------
run_sjaracne "$SIG_CLEAN" "output_sig_sjaracne_Driver_Inference_out_.final"
run_sjaracne "$TF_CLEAN"  "output_tf_sjaracne_Driver_Inference_out_.final"
run_sjaracne "$LNC_CLEAN" "output_nc_sjaracne_Driver_Inference_out_.final"

echo "All SJARACNe runs completed successfully!"
