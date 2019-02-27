
while IFS= read -r f; do
	len=${#f}
	soxi -D $f >> list_file_duration.txt
done < "$1"
