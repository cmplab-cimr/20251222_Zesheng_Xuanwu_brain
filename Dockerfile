FROM continuumio/anaconda3:latest
LABEL authors.name="Shiyu Liu" authors.email="shiyu.liu@cimrbj.ac.cn"

RUN conda update conda && \
    conda config --add channels conda-forge && \
    conda install -c conda-forge tqdm && \
    conda clean -afy
RUN pip install tensorflow tensorflow-probability xlsxwriter && \
    pip cache purge

#ENTRYPOINT ["python"]
ENTRYPOINT ["bash"]
