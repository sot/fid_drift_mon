# Licensed under a 3-clause BSD style license - see LICENSE.rst
from setuptools import setup

entry_points = {
    "console_scripts": [
        "fid_drift_mon_update_fid_stats = fid_drift_mon.update_fid_stats:main",
        "fid_drift_mon_plot_drift = fid_drift_mon.plot_drift:main",
        "fid_drift_mon_plot_drift_model = fid_drift_mon.plot_drift_model:main",
        "fid_drift_mon_plot_command_vs_telem = fid_drift_mon.plot_command_vs_telem:main",
    ]
}

setup(
    name="fid_drift_mon",
    description="Fid light drift monitor",
    author="Tom Aldcroft",
    author_email="taldcroft@cfa.harvard.edu",
    packages=["fid_drift_mon"],
    tests_require=["pytest"],
    use_scm_version=True,
    setup_requires=["setuptools_scm", "setuptools_scm_git_archive"],
    entry_points=entry_points,
    install_package_data=True,
    package_data={"fid_drift_mon": ["data/index_template.html", "task_schedule.cfg"]},
)
