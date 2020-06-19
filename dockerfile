FROM python:3.8

RUN apt-get update && apt-get install -y r-base cmake
RUN R -e "install.packages('gam', dependencies=TRUE, repos='http://cran.rstudio.com/')"

RUN pip install --trusted-host pypi.python.org cmake
RUN pip install --trusted-host pypi.python.org dash==1.7.0
RUN pip install --trusted-host pypi.python.org dash_html_components==1.0.2
RUN pip install --trusted-host pypi.python.org dash-bootstrap-components==0.8.2 
RUN pip install --trusted-host pypi.python.org plotly==4.3.0
RUN pip install --trusted-host pypi.python.org flask==1.1.1
RUN pip install --trusted-host pypi.python.org flask_caching==1.8.0
RUN pip install --trusted-host pypi.python.org gunicorn==19.9.0

RUN pip install --trusted-host pypi.python.org numba==0.48.0
RUN pip install --trusted-host pypi.python.org umap-learn==0.3.0
RUN pip install --trusted-host pypi.python.org leidenalg==0.7.0
RUN pip install --trusted-host pypi.python.org git+https://github.com/theislab/anndata.git
RUN pip install --trusted-host pypi.python.org scanpy==1.4.5.post3
RUN pip install --trusted-host pypi.python.org zarr==2.4.0
RUN pip install --trusted-host pypi.python.org git+https://github.com/dpeerlab/Palantir.git@6168f60c0ea416d87d622fccb0c0709dde2cc88c
RUN pip install --trusted-host pypi.python.org bbknn==1.3.7

RUN pip install --trusted-host pypi.python.org shortuuid==0.5.0 
RUN pip install --trusted-host pypi.python.org filelock==3.0.12 
RUN pip install --trusted-host pypi.python.org colour 


WORKDIR . /MiCV
ADD ./src/index.py ./
ADD ./src/app.py ./
ADD ./src/layouts.py ./
ADD ./src/helper_functions.py ./
ADD ./src/annotation ./annotation
ADD ./src/processing ./processing
ADD ./src/assets ./assets
ADD ./src/exporting ./exporting
ADD ./src/markergenes ./markergenes
ADD ./src/plotting ./plotting
ADD ./src/pseudotime ./pseudotime
ADD ./src/inputoutput ./inputoutput


RUN mkdir -p /run/MiCV
RUN mkdir -p /srv/www/MiCV/cache/
RUN mkdir -p /srv/www/MiCV/selected_datasets/

ADD ./gene_information_tables/*.csv /srv/www/MiCV/cache/
ADD ./selected_datasets/*.h5ad /srv/www/MiCV/selected_datasets/
ADD ./src/assets/ /srv/www/MiCV/assets/

EXPOSE 8000
EXPOSE 8050

CMD ["gunicorn", "-b", "0.0.0.0:8000", "--workers=4", "--max-requests=5", "--pid", "/run/MiCV/MiCV.pid", "index:server"]
#CMD ["python", "index.py"]