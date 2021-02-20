from setuptools import setup
from version import get_git_version

setup(name='vlass_tools',
      version=get_git_version(),
      url='http://github.com/caseyjlaw/vlass',
      packages=['vlass'],
      requirements=['astropy', 'numpy'],
      zip_safe=False)
