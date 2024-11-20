if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input.fasta> <output.txt>"
    exit 1
fi
fasta_file=$1
output_txt=$2
awk -v OFS="\t" 'BEGIN {print "GeneID\tGeneLength"}
    /^>/ { if (NR > 1) print prev_id "\t" prev_len; prev_id=$1; prev_len=0; next }
    {
        prev_len += length($0)
    }
    END { print prev_id "\t" prev_len }' $fasta_file > $output_txt
sed -i 's/>//g' "$output_txt"
echo "TXT file has been created: $output_txt"
