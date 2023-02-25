cd /home/hychimaa/workspace/molpc_imp/result/all/

for file in *
do
	echo ${file%.*}
	/home/hychimaa/workspace/molpc_imp/src/analysis/MMalign /home/hychimaa/workspace/molpc_imp/data/experimental/$file $file -outfmt 2 > ../tmscores/${file%.*}'.txt'
done
