from setuptools import find_packages, setup
from distutils.command.sdist import sdist

with open('README') as f:
  readme = f.read()

with open('LICENSE') as f:
  license = f.read()


setup(
  name='uzkChemTem',
  version='0.0.3',
  description='Python package to read TEM images using the bioformats package',
  url='https://github.com/DomiDre/uzkChemTem',
  author='Dominique Dresen',
  author_email='dominique.dresen@uni-koeln.de',
  license=license,
  long_description=readme,
  python_requires='>2.7, !=3.0.*, !=3.1.*, !=3.2.*, !=3.3.*',
  platforms=['Linux'],
  package_dir={'uzkChemTem': 'uzkChemTem'},
  packages=find_packages(
    exclude=(
      'tests',
      'examples'
      )
  ),
  include_package_data=True,
  cmdclass={'sdist': sdist},
  zip_safe=True,
  keywords='tem bioformats uzk chemistry'
)
