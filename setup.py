from setuptools import setup, find_packages

setup(
    name='CMDB_CNV',  # This is the name of your package
    version='0.1.0',  # The version of your package (you can update this as you release)
    description='CNV data normalization and visualization tools',  # Short description
    author='Your Name',  # Replace with your name or organization
    author_email='you@example.com',  # Replace with your email
    packages=find_packages(),  # This automatically finds your CMDB_CNV/CMDB_CNV package and includes it
    install_requires=[
        'anndata',  # List your package dependencies here
        'numpy',    # Add other dependencies as needed
        'scipy',    # Add any other required libraries
    ],
    classifiers=[  # Optional: This helps people find your package based on certain criteria
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.7',  # Minimum Python version required
)
