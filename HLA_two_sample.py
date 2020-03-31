import subprocess

# Using these paths for convenience
# phlat-1.0
HLA = 'python PHLAT.py'
HLA_INDEX = '~/shared/index4phlat'
HLA_PATH = "~/phlat"
BOWTIE2 = 'bowtie2'

# This function will predict HLAs using PHLAT and the two samples
def HLA_pipeline(sample1_normal, sample2_normal, sample1_cancer, sample2_cancer, threads):
    cmd1 = "{} -1 {} -2 {} -index {} -b2url {} -tag normal -e {} -o . -p {}".format(HLA,
                                                                                    sample1_normal,
                                                                                    sample2_normal,
                                                                                    BOWTIE2,
                                                                                    HLA_PATH,
                                                                                    threads)
    cmd2 = "{} -1 {} -2 {} -index {} -b2url {} -tag cancer -e {} -o . -p {}".format(HLA,
                                                                                    sample1_cancer,
                                                                                    sample2_cancer,
                                                                                    BOWTIE2,
                                                                                    HLA_PATH,
                                                                                    threads)
    p1 = subprocess.Popen(cmd1, shell=True)
    p2 = subprocess.Popen(cmd2, shell=True)
    p1.wait()
    p2.wait()
