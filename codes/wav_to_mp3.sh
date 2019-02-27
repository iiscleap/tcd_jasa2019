
while IFS= read -r f; do
	echo $f
	len=${#f}
	ffmpeg -i $f -codec:a libmp3lame -b:a 320k ${f:0:len-4}.mp3 </dev/null
	#ffmpeg -i $f -codec:a libfdk_aac -profile:a aac_he -b:a 320k ${f:0:len-4}.m4a
	#lame -q0 -b128 $f ${f:0:len-4}.mp3
done < "$1"

