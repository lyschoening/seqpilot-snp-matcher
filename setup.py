# encoding: utf-8
from distutils.core import setup

with open('requirements.txt', 'r') as f:
    requirements = f.read().split("\n")

setup(
    name='seqpilot-snp-matcher',
    version='1.0.0',
    packages=['seqpilot_snp_matcher'],
    url='',
    license='free for all',
    author=u'Lars Schöning',
    author_email='lars@lyschoening.de',
    description='PDF reports for comparing identifying SNPs in tables exported from JSI Sequence Pilot',
    entry_points={
        'console_scripts': [
            'seqpilot-snp-matcher = seqpilot_snp_matcher.report:main',
        ]
    },
    install_requires=requirements
)

