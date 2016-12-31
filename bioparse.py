from Bio.Seq import Seq


class Sequence(object):

    def __init__(self):
        self.seq_str = ""

    def get_seq(self, string_value):
        my_seq = Seq(string_value)
        return my_seq

