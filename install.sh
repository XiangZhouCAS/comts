chmod +x comts
chmod +x ./bins/com
chmod +x ./bins/custom
comts_path=$(dirname "$(readlink -f "$0")")
echo 'export PATH="$PATH:'"$comts_path"'"' >> ~/.bashrc
bins_dir=$comts_path/$(echo bins)
sed -i "1s|.*|R_SCRIPTS_DIR=$bins_dir|" $comts_path/comts
options(repos=structure(c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")))
if (!require(optparse,quietly = TRUE)) {
  install.packages("optparse")
  library(optparse)
} else {
  library(optparse)}
  if (!require(magrittr)) {
  install.packages("magrittr")
  library(magrittr)
} else {
  library(magrittr)}
if (!require(dplyr)) {
  install.packages("dplyr")
  library(dplyr)
} else {
  library(dplyr)}
if (!require(tidyr)) {
  install.packages("tidyr")
  library(tidyr)
} else {
  library(tidyr)}
if (!require(data.table)) {
  install.packages("data.table")
  library(data.table)
} else {
  library(data.table)}
if(!require(ggplot2)){
  install.packages("ggplot2")
  library(ggplt2)
}else{library(ggplot2)}
