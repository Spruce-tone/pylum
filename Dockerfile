# FROM continuumio/miniconda3:latest

# RUN apt-get update -y; apt-get upgrade -y
# RUN apt-get update -y; apt-get upgrade -y; apt-get install -y vim-tiny vim-athena ssh

# COPY environment.yml environment.yml

# RUN conda env create -f environment.yml
# RUN echo "alias l='ls -lah'" >> ~/.bashrc
# RUN echo "source activate coin" >> ~/.bashrc

# ENV CONDA_EXE /opt/conda/bin/conda
# ENV CONDA_PREFIX /opt/conda/envs/coin
# ENV CONDA_PYTHON_EXE /opt/conda/bin/python
# ENV CONDA_PROMPT_MODIFIER (coin)
# ENV CONDA_DEFAULT_ENV coin
# ENV PATH /opt/conda/envs/coin/bin:/opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

# FROM continuumio/miniconda3:latest
# RUN apt-get update -y; apt-get upgrade -y
# WORKDIR /spruce
# COPY environment.yml environment.yml
# RUN conda env create -f environment.yml

# WORKDIR /root/.jupyter
# RUN echo jupyter notebook --generate-config

# # jupyter server setting
# RUN echo '\
# from IPython.lib import passwd \n\
# password = passwd("dydwo789") \n\
# c.NotebookApp.password = password' >> /root/.jupyter/jupyter_notebook_config.py
# RUN echo "c.NotebookApp.ip = '0.0.0.0'" >> /root/.jupyter/jupyter_notebook_config.py
# RUN echo "c.NotebookApp.password_require=True" >> /root/.jupyter/jupyter_notebook_config.py
# RUN echo "c.NotebookApp.allow_root = True" >> /root/.jupyter/jupyter_notebook_config.py
# RUN echo "c.NotebookApp.open_browser = False" >> /root/.jupyter/jupyter_notebook_config.py
# RUN echo "c.NotebookApp.port=8888" >> /root/.jupyter/jupyter_notebook_config.py

# WORKDIR /spruce
# # RUN echo "conda activate pylum" >> ~/.bashrc
# RUN echo "conda activate pylum"
# CMD ["jupyter lab","--ip=0.0.0.0","--no-browser","--allow-root"]
# # # ENTRYPOINT ["jupyter", "lab","--ip=0.0.0.0","--allow-root"]

# FROM ubuntu:20.04
FROM ubuntu:latest
LABEL maintainer='spruce'
RUN apt-get update -y; apt-get upgrade -y; apt-get install -y wget \
    && wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && bash Miniconda3-latest-Linux-x86_64.sh -b -p /miniconda \
    && rm Miniconda3-latest-Linux-x86_64.sh

ENV PATH=$PATH:/miniconda/condabin:/miniconda/bin
ENV JUPYTER_ENABLE_LAB yes
ENV ENV_NAME pylum

WORKDIR /spruce
COPY environment.yml environment.yml
RUN conda env create -f environment.yml

SHELL ["sh", "-c"]
RUN conda run -n $ENV_NAME python -m ipykernel install --name $ENV_NAME --display-name [$ENV_NAME::py38]

SHELL ["/bin/bash", "-c"]
RUN conda update -y conda\ 
    && conda init \
    && echo 'conda activate $ENV_NAME' >> ~/.bashrc \
    && echo 'jupyter lab --allow-root --no-browser --ip=0.0.0.0 --notebook-dir=/spruce'>> ~/.bashrc