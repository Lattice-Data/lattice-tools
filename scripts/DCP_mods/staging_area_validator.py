"""
Runs a pre-check of a metadata files to identify issues that 
might cause the snapshot or indexing processes to fail.
Adapted from https://github.com/DataBiosphere/hca-import-validation/blob/main/hca/staging_area_validator.py
"""
import json
import logging
import os
import requests
import sys
import uuid
from typing import (
    MutableMapping,
    TypeVar,
)
from jsonschema import (
    Draft7Validator,
    FormatChecker,
    ValidationError,
)


T = TypeVar('T')
JSON = MutableMapping[str, T]
date_format = '%Y-%m-%dT%H:%M:%S.%fZ'

schemas = {}


def validate_files(path: str, reporting, output_dir) -> None:
    directory = f'{output_dir}/{reporting["dataset"]}'
    validate_file_fn = globals()[f'validate_{path}_file']
    for r, d, f in os.walk(os.path.join(directory, path)):
        for file_name in f:
            full_path = os.path.join(r, file_name)
            validate_file_fn(full_path, reporting)


def validate_links_file(file_name, reporting) -> None:
    type_file = '/'.join(file_name.split('/')[1:])
    # Expected syntax: links/{bundle_uuid}_{version}_{project_uuid}.json
    if file_name.count('_') != 2:
        reporting['file_errors'].add(type_file)
        logging.error(f"{type_file} has {file_name.count('_')} _'s in name, expecting 2")
    if not file_name.endswith('.json'):
        reporting['file_errors'].add(type_file)
        logging.error(f'{type_file} does not end in .json')
    _, _, project_uuid = file_name[:-5].split('_')
    file_json = json.load(open(file_name))
    validate_file_json(file_json, file_name, reporting)
    for link in file_json['links']:
        link_type = link['link_type']
        if link_type == 'process_link':
            add_metadata_file(entity_id=link['process_id'],
                                   entity_type=link['process_type'],
                                   project_uuid=project_uuid,
                                   category='process',
                                   reporting=reporting)
            for category in ('input', 'output', 'protocol'):
                for entity in link[f'{category}s']:
                    if entity.get(f'{category}_type'):
                        entity_type = entity[f'{category}_type']
                    else:
                        reporting['file_errors'].add(type_file)
                        logging.error(f'{type_file}, {category} missing {category}_type')
                        entity_type = 'missing-type'
                    if entity.get(f'{category}_id'):
                        entity_id = entity[f'{category}_id']
                    else:
                        reporting['file_errors'].add(type_file)
                        logging.error(f'{type_file}, {category} missing {category}_id')
                        entity_id = 'missing-id'
                    add_metadata_file(entity_id=entity_id,
                                           entity_type=entity_type,
                                           project_uuid=project_uuid,
                                           category=category,
                                           reporting=reporting)
        elif link_type == 'supplementary_file_link':
            if link['entity']['entity_type'] != 'project':
                reporting['file_errors'].add(type_file)
                logging.error(f'{type_file}, supplementary_file not linked to entity_type project')
            if link['entity']['entity_id'] != project_uuid:
                reporting['file_errors'].add(type_file)
                logging.error(f'{type_file}, supplementary_file not linked to project entity_id {project_uuid}')
            for entity in link['files']:
                entity_type = entity['file_type']
                entity_id = entity['file_id']
                add_metadata_file(entity_id=entity_id,
                                   entity_type=entity_type,
                                   project_uuid=project_uuid,
                                   category='supplementary',
                                   reporting=reporting)
    if project_uuid not in reporting['metadata_files']:
        add_metadata_file(entity_id=project_uuid,
                           entity_type='project',
                           project_uuid=project_uuid,
                           category='project',
                           reporting=reporting)


def add_metadata_file(entity_id: str,
                      entity_type: str,
                      project_uuid: str,
                      category: str,
                      reporting
                      ) -> None:
    try:
        file = reporting['metadata_files'][entity_id]
    except KeyError:
        reporting['metadata_files'][entity_id] = {
            'name': None,
            'entity_id': entity_id,
            'entity_type': entity_type,
            'version': None,
            'project': {project_uuid},
            'category': {category},
            'found_metadata': False,
            'found_descriptor': False,
        }
    else:
        file['project'].add(project_uuid)
        file['category'].add(category)


def validate_metadata_file(file_name, reporting) -> None:
    type_file = '/'.join(file_name.split('/')[1:])
    # Expected syntax: metadata/{metadata_type}/{metadata_id}_{version}.json
    metadata_type, metadata_file = file_name.split('/')[-2:]
    if metadata_file.count('_') != 1:
        reporting['file_errors'].add(type_file)
        logging.error(f"{type_file} has {metadata_file.count('_')} _'s in name, expecting 1")
    if not metadata_file.endswith('.json'):
        reporting['file_errors'].add(type_file)
        logging.error(f'{type_file} does not end in .json')
    metadata_id, metadata_version = metadata_file.split('_')
    file_json = json.load(open(file_name))
    validate_file_json(file_json, file_name, reporting)
    if file_json.get('provenance'):
        if metadata_id != file_json['provenance']['document_id']:
            reporting['file_errors'].add(type_file)
            logging.error(f"{type_file} id in file name is {metadata_id} but does not match provenance.document_id {file_json['provenance']['document_id']}")
    if reporting['metadata_files'].get(metadata_id):
        metadata_file = reporting['metadata_files'][metadata_id]
        metadata_file['name'] = file_name
        metadata_file['version'] = metadata_version
        metadata_file['found_metadata'] = True
        if metadata_type.endswith('_file'):
            metadata_file['data_file_name'] = file_json['file_core']['file_name']
    else:
        reporting['extra_files'].append(type_file)


def validate_descriptors_file(file_name, reporting) -> None:
    type_file = '/'.join(file_name.split('/')[1:])
    # Expected syntax: descriptors/{metadata_type}/{metadata_id}_{version}.json
    metadata_type, metadata_file = file_name.split('/')[-2:]
    if metadata_file.count('_') != 1:
        reporting['file_errors'].add(metadata_file)
        logging.error(f"{metadata_file} has {metadata_file.count('_')} _'s in name, expecting 1")
    if not metadata_file.endswith('.json'):
        reporting['file_errors'].add(metadata_file)
        logging.error(f'{metadata_file} does not end in .json')
    metadata_id, metadata_version = metadata_file.split('_')
    file_json = json.load(open(file_name))
    validate_file_json(file_json, file_name, reporting)
    file_name = file_json['file_name']
    reporting['names_to_id'][file_name] = metadata_id
    if reporting['metadata_files'].get(metadata_id):
        metadata_file = reporting['metadata_files'][metadata_id]
        metadata_file['found_descriptor'] = True
        metadata_file['crc32c'] = file_json['crc32c']
        version = metadata_file['version']
        if version != metadata_version:
            reporting['file_errors'].add(metadata_file)
            logging.error(f'{metadata_file} version in file name is {metadata_version} but does not match corresponding metadata file version {version}')
    else:
        reporting['extra_files'].append(type_file)


def validate_file_json(file_json: JSON, file_name: str, reporting) -> None:
    type_file = '/'.join(file_name.split('/')[1:])
    if file_json.get('describedBy'):
        if file_json['describedBy'] in schemas:
            schema = schemas[file_json['describedBy']]
        else:
            try:
                schema = download_schema(file_json['describedBy'])
                schemas[file_json['describedBy']] = schema
            except json.decoder.JSONDecodeError as e:
                reporting['file_errors'].add(type_file)
                logging.error(f"{type_file}, failed to parse schema JSON from {file_json['describedBy']}")
    else:
        reporting['file_errors'].add(type_file)
        logging.error(f'{type_file}, describedBy is required')
    if 'schema' in locals():
        validator = Draft7Validator(schema)
        errors = validator.iter_errors(file_json)
        for e in errors:
            reporting['file_errors'].add(type_file)
            logging.error(f'{type_file}, {e.message}')


def download_schema(schema_url: str) -> JSON:
    response = requests.get(schema_url, allow_redirects=False)
    response.raise_for_status()
    return response.json()


def check_results(reporting):
    for metadata_id, metadata_file in reporting['metadata_files'].items():
        if not metadata_file['found_metadata']:
            reporting['file_errors'].add(metadata_file['entity_id'])
            logging.error(f"{metadata_file['entity_type']} {metadata_file['entity_id']} metadata file not found")
        if metadata_file['entity_type'].endswith('_file') and not metadata_file['found_descriptor']:
            reporting['file_errors'].add(metadata_file['entity_id'])
            logging.error(f"{metadata_file['entity_type']} {metadata_file['entity_id']} descriptor file not found")


#def write_error(message, output_dir):
#    with open(f'{output_dir}/validation_errors.txt', 'a') as outfile:
#        outfile.write(message + '\n')


def dcp_validation(dataset, output_dir):
    #logging.basicConfig(filename=f'{output_dir}/validation.log', level=logging.INFO)
    logging.info('validating ' + dataset)
    reporting = {
        'dataset': dataset,
        'names_to_id': {},
        'metadata_files': {},
        'file_errors': set(),
        'extra_files': [],
    }
    for t in ['links','metadata','descriptors']:
        logging.info('validating ' + t)
        validate_files(t, reporting, output_dir)
    check_results(reporting)

    for file_name in reporting['extra_files']:
        reporting['file_errors'].add(file_name)
        logging.error(f'{file_name}, not a part of a subgraph in links files')

    if not reporting['file_errors']:
        logging.info('No DCP schema validation errors')

    return(len(reporting['file_errors']))
