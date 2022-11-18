#pwd
#readarray -t a < /Users/cfgp/Sun-Lab-CHiP-seq/random.txt
IFS=$'\n' read -d '' -r -a lines < /Users/cfgp/Sun-Lab-CHiP-seq/random.txt
#printf '%s\n' "${lines[@]}"
IFS=$'\n' read -d '' -r -a anotherlist < /Users/cfgp/Sun-Lab-CHiP-seq/another.txt

VAR=""
placeholder="~directory/to/tagdirs/"
for ELEMENT in ${lines[@]}; do
  VAR+="${placeholder}"
  VAR+="${ELEMENT} "
done

for line in ${anotherlist[@]}; do
    VAR+="${placeholder}"
    VAR+="${line} "
done
echo "$VAR"

#printf '%s\n' "${lines[@]}"