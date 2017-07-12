
for entry in `ls $search_dir`; do
    echo `file $entry`
done