# This Python file uses the following encoding: utf-8
import json

class LeeDistance():
    def __init__(self,file_name):
        with open(file_name, 'r') as f:
            self.n_list = json.load(f)

    def get_min_lee_dist_for_k1_k2(self,n,k1,k2):
        lee_dict = self.n_list[str(n)]
        key = str(k1)+"-"+str(k2)
        return lee_dict[key]


# if __name__ == "__main__":
#     pass