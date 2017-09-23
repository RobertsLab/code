#!/bin/bash

# This script was written to address a specific problem that Rhonda was having
# where multiple GO terms are present in a single field, but she needs them
# on separate lines.
# See this GitHub issue for explanation of the problem:
# https://github.com/sr320/LabDocs/issues/654

# Takes a tab-delimited input file that contains multiple GO terms, separated by a ";<space>", in a single field.

# GO terms in input file MUST BE IN LAST FIELD (i.e. column).

# Before executing this script, the following need to be changed (look for them below):

# /path/to/input/file
# /path/to/new/filename
# /path/to/output_file


# Set variables for files
# input_file is the initial, "problem" file
# tmp_file is an intermediate file that most of the program works upon
# output_file is the final file produced by the script
input_file="/path/to/input/file"
tmp_file="/path/to/new/filename"
output_file="/paht/to/output/file"

# sed command substitutes the "; " sequence to a tab and writes the new format to a new file.
# This character sequence is how the GO terms are delimited in their field.
sed $'s/; /\t/g' "$input_file" > "$tmp_file"

# Identify first field containing a GO term.
# Search file with grep for "GO:" and pipe to awk.
# Awk sets tab as field delimiter (-F'\t'), runs a for loop that looks for "GO:" (~/GO:/), and then prints the field number).
# Awk results are piped to sort, which sorts unique by number (-ug).
# Sort results are piped to head to retrieve the lowest value (i.e. the top of the list; "-n1").
begin_goterms=$(grep "GO:" "$tmp_file" | awk -F'\t' '{for (i=1;i<=NF;i++) if($i ~/GO:/) print i}' | sort -ug | head -n1)

# Save value of field number that indicates last field before GO terms.
end_fixed_fields=$(($begin_goterms-1))

# While loop to process each line of the input file.
# Initial while loop statement is written to handle the last line of the file, even if it doesn't end with a newline character.
while IFS='' read -r line || [ -n "$line" ]
	do
	
	# Send contents of the current line to awk.
	# Set the field separator as a tab (-F'\t') and print the number of fields in that line.
	# Save the results of the echo/awk pipe (i.e. number of fields) to the variable "max_field".
	max_field=$(echo "$line" | awk -F'\t' '{print NF}')

	# Send contents of current line to cut.
	# Cut all fields (i.e. retain those fields) up to (not including) the first field with with a GO term.
	# Save the results of the echo/cut pipe (i.e. fields 1-12) to the variable "fixed_fields"
	fixed_fields=$(echo "$line" | cut -f1-$end_fixed_fields)

	# Since not all the lines contain the same number of fields (e.g. may not have GO terms),
	# evaluate the number of fields in each line to determine how to handle current line.

	# If the value in max_field is less than the field number where the GO terms begin,
	# then just print the current line (%s) followed by a tab (\t) and newline (\n).
	if (( "$max_field" < "$begin_goterms" ))
		then printf "%s\t\n" "$line"
			else

			# Send contents of current line (which contains GO terms) to cut.
			# Cut fields (i.e. retain those fields) to whatever the last field is in the curent line.
			# Save the results of the echo/cut pipe (i.e. all the GO terms fields) to the variable "goterms".
			goterms=$(echo "$line" | cut -f"$begin_goterms"-"$max_field")
			
			# Assign values in the variable "goterms" to a new indexed array (called "array"), 
			# with tab delimiter (IFS=$'\t')
			IFS=$'\t' read -r -a array <<<"$goterms"
			
			# Iterate through each element of the array.
			# Print the first 12 fields (i.e. the fields stored in "fixed_fields") followed by a tab (%s\t).
			# Print the current element in the array (i.e. the current GO term) followed by a new line (%s\n).
			for element in "${!array[@]}"	
				do printf "%s\t%s\t\n" "$fixed_fields" "${array[$element]}"
			done
	fi

# Send the input file into the while loop and send the output to a file named "rhonda_fixed.txt".
done < "$tmp_file" > "$output_file"
