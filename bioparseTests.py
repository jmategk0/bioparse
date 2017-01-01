import bioparse
import time
# import json

# from pprint import pprint
parser = bioparse.Genbank()


start_time = time.time()
# d = parser.genbank_to_dictionary(genbank_file="test_files/biopython/ls_orchid.gbk")
d = parser.genbank_to_dictionary(genbank_file="test_files/biopython/ls_orchid.gbk")
end_time = time.time()
# e = json.dumps(d)

# pprint(d)
print(type(d))
print("--- %s seconds ---" % (end_time - start_time))
# parser.genbank_to_dictionary
# --- 0.03371381759643555 seconds ---
# --- 0.03315591812133789 seconds ---
# --- 0.03240847587585449 seconds ---


# parser.genbank_to_dictionary_lite
# --- 0.03237748146057129 seconds ---
# --- 0.03287935256958008 seconds ---
# --- 0.03213834762573242 seconds ---

# diff is ~0.001 seconds. Will this matter at scale? What is more readable?
