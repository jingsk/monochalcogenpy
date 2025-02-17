from distutils.core import setup
setup(name='monochalcogen_analysis',
      version='1.0',
      #package_dir={'becqsdr': 'utils'},
      packages=['monochalcogenpy'],
      install_requires = [
        'ase>=3.23',
        'matscipy>=1.1.0',
    ]
      )
