import json
import os
import sys
import requests
from pathlib import Path
from requests.auth import HTTPBasicAuth


def get_path(search_term: str) -> os.PathLike | str:
    """
    Find path of local repos and API keys regardless of source machine.
    Returns Path when found, otherwise str "Path not found"
    """
    local_path = Path()
    
    likely_locations = [
        local_path.home() / "Documents" / "keys",
        local_path.home() / "keys",
    ]

    for place in likely_locations:
        if place.exists():
            for item in place.iterdir():
                if search_term in item.name:
                    return item

    return "Path not found"


class JIRA_API:
    """JIRA API access configuration and helper methods"""
    
    # Class variables to store config
    project_id = None
    email = None
    api_token = None
    jira_url = None
    api_url = None
    auth = None
    headers = None
    
    @classmethod
    def config(cls, config_file: str = "jira.json"):
        """
        Configure JIRA API access by loading credentials from config file.
        
        :param config_file: name of the JSON config file (default: jira.json)
        :return: None
        """
        # Find config file
        config_path = get_path(config_file)
        
        if not isinstance(config_path, Path):
            print(f"Config file '{config_file}' not found in Documents/keys/ or other standard locations")
            sys.exit(1)
        
        # Load config
        try:
            with open(config_path, 'r') as f:
                config = json.load(f)
        except json.JSONDecodeError as e:
            print(f"ERROR: Invalid JSON in config file: {e}")
            sys.exit(1)
        
        # Validate required fields
        required_fields = ['project_id', 'email', 'api_token']
        missing_fields = [field for field in required_fields if not config.get(field)]
        
        if missing_fields:
            print(f"ERROR: Missing required fields in config: {', '.join(missing_fields)}")
            sys.exit(1)
        
        # Set class variables
        cls.project_id = config['project_id']
        cls.email = config['email']
        cls.api_token = config['api_token']
        cls.jira_url = f"https://{cls.project_id}.atlassian.net"
        cls.api_url = f"{cls.jira_url}/rest/api/3/search/jql"
        cls.auth = HTTPBasicAuth(cls.email, cls.api_token)
        cls.headers = {"Accept": "application/json"}
        
        print(f"JIRA API configured for: {cls.jira_url}")
    
    @classmethod
    def get_all_issues(cls, jql: str, fields: list = None, max_results: int = 100):
        """
        Fetch all issues matching a JQL query with automatic pagination.
        
        :param jql: JQL query string (e.g., 'project = CXG')
        :param fields: List of fields to retrieve (default: all fields)
        :param max_results: Number of results per page (default: 100)
        :return: List of all issues
        """
        if not cls.api_url:
            raise RuntimeError("JIRA API not configured. Call JIRA_API.config() first.")
        
        all_issues = []
        next_page_token = None
        
        # Default fields if not specified
        if fields is None:
            fields = ['key', 'status', 'customfield_10041', 'labels']
        
        while True:
            query = {
                'jql': jql,
                'maxResults': max_results,
                'fields': ','.join(fields) if isinstance(fields, list) else fields
            }
            
            # Add nextPageToken if we have one
            if next_page_token:
                query['nextPageToken'] = next_page_token
            
            response = requests.get(
                cls.api_url,
                headers=cls.headers,
                params=query,
                auth=cls.auth
            )
            
            response.raise_for_status()  # Raise exception for bad status codes
            
            data = response.json()
            issues = data.get('issues', [])
            all_issues.extend(issues)
            
            # Check for next page token
            next_page_token = data.get('nextPageToken')
            if not next_page_token:
                break

        print(f'Found {len(all_issues)} JIRA issue')

        return all_issues

    @classmethod
    def search_issues(cls, jql: str, fields: list = None, max_results: int = 100):
        """
        Alias for get_all_issues for consistency with other APIs.
        """
        return cls.get_all_issues(jql, fields, max_results)

    @classmethod
    def get_issue(cls, issue_key: str, fields: list = None):
        """
        Fetch a single JIRA issue by its key.
        
        :param issue_key: The issue key (e.g., 'CXG-123')
        :param fields: List of fields to retrieve (default: all fields)
        :return: Dictionary containing the issue data, or None if not found
        """
        if not cls.jira_url:
            raise RuntimeError("JIRA API not configured. Call JIRA_API.config() first.")
        
        # Build the URL for single issue endpoint
        issue_url = f"{cls.jira_url}/rest/api/3/issue/{issue_key}"
        
        # Build query parameters
        params = {}
        if fields:
            params['fields'] = ','.join(fields) if isinstance(fields, list) else fields
        
        try:
            response = requests.get(
                issue_url,
                headers=cls.headers,
                params=params,
                auth=cls.auth
            )
            
            response.raise_for_status()
            
            issue_data = response.json()
            return issue_data
            
        except requests.exceptions.HTTPError as e:
            if e.response.status_code == 404:
                print(f"Issue '{issue_key}' not found")
                return None
            else:
                print(f"Error fetching issue '{issue_key}': {e}")
                raise
        except requests.exceptions.RequestException as e:
            print(f"Request error: {e}")
            raise
