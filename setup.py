from setuptools import setup

setup(
    name = 'HaSAPPy',
    packages = ['HaSAPPy'],
    version = '1.0.0',
    description = 'analysis of FW genetic screenings in Haploid cells',
    author = 'Giulio Di Minin',
    author_email = 'giulio@diminin.it',
    url = 'https://github.com/gdiminin/HaSAPPy',
    license = 'MIT',
    download_url = 'https://github.com/gdiminin/HaSAPPy.git',
    keywords = [],
    classifiers = [],
    install_requires=['numpy',
      'HTSeq',
      'matplotlib',
      'pandas',
      'scipy',
      'xlsxwriter',
      'sklearn'],
    include_package_data = True,
    )
