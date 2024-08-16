from setuptools import setup, find_packages

setup(
    name='parse_matlab_code',
    version='0.0.1',
    author='redpigkiller',
    author_email='hongrunlu@outlook.com',
    description='Parse MATLAB code',
    url='https://github.com/redpigkiller/parse_matlab_code',
    packages=find_packages(),
    package_data={
        'parse_matlab_code': ['grammar/*.yaml', 'grammar/*.lark'],
    },
    include_package_data=True,
)