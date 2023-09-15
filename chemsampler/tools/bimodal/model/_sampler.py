import os
from sample import Sampler
import sys
import csv

root = os.path.dirname(os.path.abspath(__file__))


class _sampler:
    def __init__(self):
        self.cwd = os.getcwd()
        self.exec_folder = os.path.join(root)
        self.temperature = 0.7
        self.fold = [1]
        self.epoch = [9]

    def sampler_ForwardRNN(self, n, data_file):
        experiment_name = "BIMODAL_random_1024"
        s = Sampler(experiment_name)
        samples = s.sample(
            N=n,
            T=self.temperature,
            fold=self.fold,
            epoch=self.epoch,
            valid=True,
            novel=True,
            unique=True,
            write_csv=False,
        )
        with open(data_file, "w") as f:
            write = csv.writer(f)
            write.writerow(["smiles"])
            for sample in samples[0]:
                write.writerow([sample])

    # def sampler_BIMODAL(self, n):
    #     experiment_name = 'BIMODAL_random_1024'
    #     s = Sampler(experiment_name)
    #     s.sample(N=n, T=self.temperature, fold=self.fold, epoch= self.epoch,
    #              valid=True, novel=True, unique=True, write_csv=True)


if __name__ == "__main__":
    n = sys.argv[1]
    data_file = sys.argv[2]
    _sampler().sampler_ForwardRNN(n, data_file)
