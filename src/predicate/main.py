from __future__ import absolute_import
import click
import json
import yaml
from yaml import Loader
from omegaconf import OmegaConf
from omegaconf.errors import ConfigKeyError
from .utils import *
from .predicate import main as predicate_main
from .config import default_config


def validate_config(config_file):
    # Check if file is json or yaml file
    if not config_file.endswith('.yml'):
        raise ValueException('Invalid config file. Only .yml files are supported.')


def validate_mandatory_fields(cfg):
    if 'metabolic_network' not in cfg:
        raise ValueError('Missing mandatory field in config: metabolic_network')
    if 'model' not in cfg.metabolic_network:
        raise ValueError('Missing mandatory field in config: metabolic_network.model')
    if 'objective' not in cfg.metabolic_network:
        raise ValueError('Missing mandatory field in config: metabolic_network.objective')
    if 'reference_files' not in cfg:
        raise ValueError('Missing mandatory field in config: reference_files')
    if 'fasta_sequences' not in cfg.reference_files:
        raise ValueError('Missing mandatory field in config: reference_files.fasta_sequences')
    if 'protein_sequences' not in cfg.reference_files:
        raise ValueError('Missing mandatory field in config: reference_files.protein_sequences')


def filter_empty(cfg):
    cfg['metabolic_network']['additional_metabolites'] = [x for x in cfg['metabolic_network']['additional_metabolites'] if len(x) > 0]
    cfg['structural_proteins'] = [x for x in cfg['structural_proteins'] if len(x) > 0]
    return cfg


@click.command()
@click.argument('config_file', type=click.Path(exists=True))
def predicate(config_file):
    validate_config(config_file)
    cfg = OmegaConf.load(config_file)
    validate_mandatory_fields(cfg)
    default = default_config()
    OmegaConf.set_struct(default, True)
    try:
        cfg = OmegaConf.merge(default, cfg)
    except ConfigKeyError as e:
        ee = 'Invalid key in config file: ' + str(e)
        raise KeyError(ee)
    cfg = filter_empty(cfg)
    print("PREDICATE Settings:")
    print("------------------------------------------")
    print(OmegaConf.to_yaml(cfg))
    print("------------------------------------------")
    predicate_main(cfg)

if __name__ == "__main__":
    predicate()
