from Bio import SeqIO


class BaseParser(object):

    def __init__(self):
        self.filename = ""
        self.data = ""

    def rename_dictionary_keys(self, dictionary, dictionary_key_map):
        """

        :param dictionary: dictionary with keys that we want to rename.
        :param dictionary_key_map: dictionary that maps old key names (keys) with new key names (value).
        :return: dictionary where keys have been replaced with the values from dictionary_key_map
        """

        for key, value in dictionary_key_map.items():
            dictionary[value] = dictionary.pop(key)
        return dictionary

    def list_of_dictionary_records_to_dictionary(self, list_of_dictionaries, primary_key='id'):
        """

        :param list_of_dictionaries: A list where each element is a dictionary.
        :param primary_key: A unique identifier key contained in each list element dictionary.
        :return: A new dictionary where each key is the primary_key and value is a dictionary from the list.
        """

        dictionary = {}

        for record in list_of_dictionaries:
            dictionary[record[primary_key]] = record
        return dictionary

    def remove_keys_from_dictionary(self, dictionary, list_of_keys):

        if isinstance(list_of_keys, str):
            list_of_keys = [list_of_keys]

        for key_name in list_of_keys:
            dictionary.pop(key_name)
        return dictionary

    def object_to_dictionary(self, value):

        return value.__dict__  # value.__dict__  # vars(value)

    def genbank_to_dictionary(self, genbank_file):

        raise NotImplementedError

    def fasta_to_dictionary(self, fasta_filename):

        raise NotImplementedError

    def fastq_to_dictionary(self, fastq_filename):

        raise NotImplementedError

    def fastq_to_fasta(self, fastq_filename, fasta_filename):

        raise NotImplementedError

    def genbank_to_fasta(self, genbank_filename, fasta_filename):

        raise NotImplementedError

    def fasta_to_genbank(self, genbank_filename, fasta_filename):

        raise NotImplementedError


class BioPython(BaseParser):

    def clean_featurelocation(self, featurelocations):

        if isinstance(featurelocations, list):
            locations = []
            for featurelocation in featurelocations:
                featurelocation_dict = self.object_to_dictionary(value=featurelocation)
                featurelocation_dict["_start"] = int(featurelocation_dict["_start"])
                featurelocation_dict["_end"] = int(featurelocation_dict["_end"])
                locations.append(featurelocation_dict)
            return_value = locations
        else:
            featurelocation_dict = self.object_to_dictionary(value=featurelocations)
            featurelocation_dict["_start"] = int(featurelocation_dict["_start"])
            featurelocation_dict["_end"] = int(featurelocation_dict["_end"])
            return_value = featurelocation_dict
        return return_value

    def clean_seqfeatures(self, seqfeatures):
        feature_list = []

        for feature in seqfeatures:
            feature_dict = self.object_to_dictionary(value=feature)

            if feature_dict['id'] == '<unknown id>':
                # Set "unknown id" to None to be more pythonic and better support json serialisation.
                feature_dict['id'] = None
                feature_dict['location'] = self.clean_featurelocation(feature_dict['location'])

            feature_list.append(feature_dict)
        return feature_list

    def clean_annotation_references(self, annotation_references):

        references = []
        for reference in annotation_references:
            reference_dict = self.object_to_dictionary(value=reference)
            if reference_dict["location"]:
                reference_dict["location"] = self.clean_featurelocation(reference_dict["location"])
            references.append(reference_dict)
        return references

    def clean_seq(self, biopython_seq):
        biopython_seq_dict = self.object_to_dictionary(value=biopython_seq)
        biopython_seq_dict["alphabet"] = biopython_seq_dict["alphabet"].letters
        return biopython_seq_dict

    def genbank_to_dictionary(self, genbank_file):

        clean_genbank_records = []

        for genbank_record_as_seq_record_obj in SeqIO.parse(genbank_file, "genbank"):

            raw_seq_record_dict = self.object_to_dictionary(value=genbank_record_as_seq_record_obj)

            raw_seq_record_dict["_seq"] = self.clean_seq(raw_seq_record_dict["_seq"])

            if raw_seq_record_dict["annotations"]["references"]:
                raw_seq_record_dict["annotations"]["references"] = self.clean_annotation_references(
                    raw_seq_record_dict["annotations"]["references"])

            if raw_seq_record_dict["features"]:
                raw_seq_record_dict["features"] = self.clean_seqfeatures(raw_seq_record_dict["features"])

            clean_genbank_records.append(raw_seq_record_dict)

        return clean_genbank_records

    def genbank_to_dictionary_lite(self, genbank_file):

        clean_genbank_records = []

        for genbank_record_as_seq_record_obj in SeqIO.parse(genbank_file, "genbank"):

            raw_seq_record_dict = genbank_record_as_seq_record_obj.__dict__

            raw_seq_record_dict["_seq"] = raw_seq_record_dict["_seq"].__dict__
            raw_seq_record_dict["_seq"]["alphabet"] = raw_seq_record_dict["_seq"]["alphabet"].letters

            if raw_seq_record_dict["annotations"]["references"]:
                references = []
                for reference in raw_seq_record_dict["annotations"]["references"]:
                    reference_dict = reference.__dict__
                    if reference_dict["location"]:
                        feature_locations = []

                        for feature_location in reference_dict["location"]:
                            feature_location_dict = feature_location.__dict__
                            feature_location_dict["_start"] = int(feature_location_dict["_start"])
                            feature_location_dict["_end"] = int(feature_location_dict["_end"])
                            feature_locations.append(feature_location_dict)
                            reference_dict["location"] = feature_locations
                    references.append(reference.__dict__)
                raw_seq_record_dict["annotations"]["references"] = references

            if raw_seq_record_dict["features"]:
                features = []
                for feature in raw_seq_record_dict["features"]:
                    feature_dict = feature.__dict__
                    if feature_dict["id"] == '<unknown id>':
                        feature_dict["id"] = None
                    feature_dict["location"] = feature_dict["location"].__dict__
                    feature_dict["location"]["_start"] = int(feature_dict["location"]["_start"])
                    feature_dict["location"]["_end"] = int(feature_dict["location"]["_end"])
                    features.append(feature_dict)
                raw_seq_record_dict["features"] = features
            clean_genbank_records.append(raw_seq_record_dict)

        return clean_genbank_records

    def fasta_to_dictionary(self, fasta_filename, raw_biopython=True):

        keys_to_remove = ["_per_letter_annotations", "annotations", "dbxrefs", "features"]
        sequences = []

        if raw_biopython:

            for seq_record_obj in SeqIO.parse(fasta_filename, "fasta"):
                sequence = self.object_to_dictionary(value=seq_record_obj)
                sequence["_seq"] = self.clean_seq(biopython_seq=sequence["_seq"])
                sequences.append(sequence)
        else:
            for seq_record_obj in SeqIO.parse(fasta_filename, "fasta"):
                sequence = self.object_to_dictionary(value=seq_record_obj)
                sequence["seq"] = str(sequence["_seq"])
                sequence.pop("_seq")
                sequence = self.remove_keys_from_dictionary(sequence, keys_to_remove)
                sequences.append(sequence)
        return sequences

    def fastq_to_dictionary(self, fastq_filename, fastq_format="fastq", raw_biopython=True):
        # Format options: ["fastq-solexa", "fastq-illumina", "fastq"]
        # SRR020192 example file

        keys_to_remove = ["annotations", "dbxrefs", "features"]
        sequences = []

        if raw_biopython:
            for seq_record_obj in SeqIO.parse(fastq_filename, fastq_format):
                sequence = self.object_to_dictionary(value=seq_record_obj)
                sequence["_seq"] = self.clean_seq(biopython_seq=sequence["_seq"])
                sequences.append(sequence)
        else:
            for seq_record_obj in SeqIO.parse(fastq_filename, fastq_format):
                sequence = self.object_to_dictionary(value=seq_record_obj)
                sequence["seq"] = str(sequence["_seq"])
                sequence.pop("_seq")
                sequence = self.remove_keys_from_dictionary(sequence, keys_to_remove)
                sequences.append(sequence)
        return sequences

    def filter_fastq(self):
        # just an example from the biopython docs 
        good_reads = (rec for rec in SeqIO.parse("SRR020192.fastq", "fastq") if min(rec.letter_annotations["phred_quality"]) >= 20)
        count = SeqIO.write(good_reads, "good_quality.fastq", "fastq")
        print("Saved %i reads" % count)

    def fastq_to_fasta(self, fastq_filename, fasta_filename):

        count = SeqIO.convert(fastq_filename, "fastq", fasta_filename, "fasta")
        return count

    def genbank_to_fasta(self, genbank_filename, fasta_filename):

        count = SeqIO.convert(genbank_filename, "genbank", fasta_filename, "fasta")
        return count

    def fasta_to_genbank(self, genbank_filename, fasta_filename):

        count = SeqIO.convert(fasta_filename, "fasta", genbank_filename, "genbank")
        return count


class BioJSON(BaseParser):

    def genbank_to_dictionary(self, genbank_file):

        raise NotImplementedError

    def fasta_to_dictionary(self, fasta_filename):

        raise NotImplementedError

    def fastq_to_dictionary(self, fastq_filename):

        raise NotImplementedError

    def fastq_to_fasta(self, fastq_filename, fasta_filename):

        raise NotImplementedError

    def genbank_to_fasta(self, genbank_filename, fasta_filename):

        raise NotImplementedError

    def fasta_to_genbank(self, genbank_filename, fasta_filename):

        raise NotImplementedError
