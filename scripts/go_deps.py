import os
import pkg_resources

REQUIRED_PACKAGES = [
    'pandas',
    'matplotlib',
    'goatools',
    'gseapy',
    'biopython'
]

for package in REQUIRED_PACKAGES:
    try:
        dist = pkg_resources.get_distribution(package)
        print('{} ({}) is installed'.format(dist.key, dist.version))
    except pkg_resources.DistributionNotFound:
        print('{} is NOT installed'.format(package))
        os.system('pip install ' + package)
