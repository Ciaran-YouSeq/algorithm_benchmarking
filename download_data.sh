#!/bin/bash

input="$1"
while IFS= read -r line
do
	wget --directory-prefix /data/reads/ "$line"
done < "$input"
