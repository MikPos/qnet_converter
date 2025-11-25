#!/usr/bin/bash


QNET_DIR="/gpfs/projects/qm_organics/posselt"

cd "$QNET_DIR"
echo "Current Directory: $(pwd)"

educts="educt.xyz.connectivity"
products="product.xyz.connectivity"
copied_educts="false"

read -p 'Type the name of Qnet directory: ' workdir

cd "$workdir"
mkdir "files_for_mod"
result_dir="$QNET_DIR/files_for_mod"

cp "relations.json" "$workdir"

# Iterate over files and directories in the working directory:
for educt_dir in "$workdir"/*_*/; do
    echo "processing: $educt_dir"
    for dir in /PATH*/; do
        echo "processing files from: $dir"
        for file in "$dir"/*.connectivity; do
            filename=$(basename "$file")
            if [ "$copied_educts" = "false" ] && [ "$filename" = "$educts" ]; then
                cp "$file" "$result_dir/$(basename "$educt_dir").xyz.connectivity"
                copied_educts="true"
            fi
            if [ "$filename" = "$products" ]; then
                cp "$file" "$result_dir/$(basename "$dir").xyz.connectivity"
                echo "$result_dir/$(basename "$dir").xyz.connectivity"
            fi
        done
    done
    copied_educts="false"
done


