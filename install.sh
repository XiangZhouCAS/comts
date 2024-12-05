chmod +x comts
chmod +x com
chmod +x custom
comts_path=$(dirname "$(readlink -f "$0")")
echo 'export PATH="$PATH:'"$comts_path"'"' >> ~/.bashrc
bins_dir=$comts_path/$(echo bins)
sed -i "1s|.*|R_SCRIPTS_DIR=$bins_dir|" $comts_path/comts
