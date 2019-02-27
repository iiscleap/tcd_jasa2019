x1='./'
x2='.wav'
x3='/Users/neeks/Desktop/Documents/work/code/matlab_codes/spkr_Change_Det/data/ibmWatson/flac/'
x4='.txt'

find $x3 -type f -maxdepth 1 -name 'LS*.flac' | while read -r line
do
        b=$(basename $line)
        len=${#b}
        echo $b
	curl -X POST -u T -u 6a8adf19-8ae5-47ea-8cc4-edf4a2d0fed2:FTNvFGZTyEEL --header "Content-Type: audio/flac" --data-binary @"$line"  "https://stream.watsonplatform.net/speech-to-text/api/v1/recognize?model=en-US_NarrowbandModel&speaker_labels=true" > "../output/${b:0:$len-5}_ibmWatson_recog.txt"
done
