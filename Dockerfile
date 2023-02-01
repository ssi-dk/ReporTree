# syntax=docker/dockerfile:1
# Some later versions of continuumio/miniconda3 cause problems (KeyError('pkgs_dirs'))
FROM continuumio/miniconda3:4.12.0
ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1
RUN conda install -c etetoolkit -c anaconda -c bioconda python=3.8 biopython=1.77 pandas=1.1.3 numpy=1.19.2 grapetree=2.1 treecluster=1.0.3 ete3 scikit-learn cgmlst-dists git --yes

COPY scripts/ /scripts/

RUN mkdir /mnt/test_data
COPY examples/Listeria/input/ /mnt/test_data/

RUN useradd -ms /bin/bash myuser
USER myuser
WORKDIR /home/myuser
RUN mkdir workdir
WORKDIR /home/myuser/workdir

# CMD ["python", "/scripts/keep_running.py"]