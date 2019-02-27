x1='./'
x2='.wav'
x3='/Users/neeks/Desktop/Documents/work/code/matlab_codes/spkr_Change_Det/data/libriSpeech/spkrChange/listTest/M_5_spkrs_8_pairs_v0/'
x4='.flac'

find $x3 -type f -maxdepth 1 -name 'LS*.wav' | while read -r line
do
        b=$(basename $line)
        len=${#b}
        echo $b
        sox $line "../flac/${b:0:len-4}$x4" rate 8000
	#printf '%s %s\n' "$line" "$x3${b:0:$len-4}$x4" >> "$create_file"
done
