from Bio import SeqIO


class Genbank(object):

    def __init__(self):
        self.filename = ""

    def object_to_dictionary(self, value):

        return value.__dict__  # value.__dict__  # vars(value)

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
