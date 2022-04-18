from setuptools import setup, find_packages

setup(
    name='FancyFlexScore',
    version='1.0',
    #packages=['FancyFlexScore'],
    url="",
    data_files = [("", ['LICENSE'])],
    license = "MIT",
    description='Flexibility score for proteins',
    author='Iria Pose, Gerard Romero Sola, Leidy Alejandra Gonzalez Molano',
    author_email='iria.pose01@estudiant.upf.edu, gerard.romero02@estudiant.upf.edu and aleja',
    install_requires = ['biopython','matplotlib','statistics','requests','pandas','seaborn', 'numpy'],
    scripts = ['__main__.py'],
    keywords = '',
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License"
        ],
     include_package_data=True
)
