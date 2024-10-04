import numpy as np
import xml.etree.ElementTree as ET
import json


# Define a function to extract desired values
def basta_lp_loader(file_path):
    # Dictionary to store the extracted values
    extracted_data = {}

    # Define the keys you want to extract with their user-friendly names
    keys_mapping = {
        "ss": "Sd",
        "cms": "Cms",
        "mms": "Mms",
        "bl": "Bl",
        "Qts": "Qts",
        "Re": "Re",
        "Rms": "Rms",
        "Qes": "Qes",
        "Le": "Le"
    }

    # Open the file and read line by line
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()

            # Stop processing when [fgprio] section is reached
            if line.startswith('[fgprio]'):
                break

            # Split each line by '=' to extract key-value pairs
            if '=' in line:
                key, value = line.split('=')
                key = key.strip()  # Clean up key
                value = value.strip()  # Clean up value

                # Check if the key is in the mapping and not already extracted
                if key in keys_mapping:
                    extracted_data[keys_mapping[key]] = float(value)

    return extracted_data


# Define a function to extract desired values from XML
def qspeaker_lp_loader(file_path):
    # Parse the XML file
    tree = ET.parse(file_path)
    root = tree.getroot()

    # Find the speaker element
    speaker = root.find('speaker')

    # Dictionary to store the extracted values
    extracted_data = {}

    # Extract the relevant attributes from the speaker element
    if speaker is not None:
        # Define the keys you want to extract
        keys = {
            "le": "Le",
            "re": "Re",
            "bl": "Bl",
            "sd": "Sd",
            "qes": "Qes",
            "qms": "Qms",
            "qts": "Qts",
            "fs": "Fs"
                
        }
        
        # Extract the values and store in the dictionary
        for xml_key, user_friendly_key in keys.items():
            extracted_data[user_friendly_key] = float(speaker.attrib.get(xml_key, None))
    
    
    Qes = extracted_data["Qes"]
    Qms = extracted_data["Qms"]
    Bl  = extracted_data["Bl"]
    Re  = extracted_data["Re"]
    Fs  = extracted_data["Fs"]
    extracted_data["Mms"] = round(Qes * Bl**2 * 1/2/np.pi/Fs/Re, 5)
    extracted_data["Cms"] = round((1/2/np.pi/Fs)**2  / 
                                  extracted_data["Mms"], 6)
    extracted_data["Rms"] = 2*np.pi*Fs*extracted_data["Mms"]/Qms
    return extracted_data



# Define a function to extract desired values from JSON
def speakerSim_lp_loader(file_path):
    # Open and load the JSON file
    with open(file_path, 'r') as file:
        data = json.load(file)

    # Define the keys you want to extract
    keys = {
        "Le": "Le",
        "Re": "Re",
        "Bl": "Bl",
        "Sd": "Sd",
        "Cms": "Cms",
        "Mms": "Mms",
        "Rms": "Rms"
    }

    # Dictionary to store the extracted values
    extracted_data = {}

    # Extract the values and store in the dictionary
    for json_key, user_friendly_key in keys.items():
        extracted_data[user_friendly_key] = float(data.get(json_key, None))

    return extracted_data


def hornResp_lp_loader(file_path):
    # Dictionary to store the extracted values
    extracted_data = {}

    # Define the keys you want to extract
    keys = {
        "Sd": "Sd",
        "Bl": "Bl",
        "Cms": "Cms",
        "Rms": "Rms",
        "Mmd": "Mms",  # Using Mms instead of Mmd for consistency
        "Le": "Le",
        "Re": "Re"
    }

    # Open the file and read line by line
    with open(file_path, 'r') as file:
        for line in file:
            # Split each line by '='
            if '=' in line:
                key, value = line.strip().split('=')
                key = key.strip()  # Clean up key
                value = value.strip()  # Clean up value

                # Check if the key is one of the desired ones and not already extracted
                if key in keys and keys[key] not in extracted_data:
                    extracted_data[keys[key]] = float(value)

    return extracted_data


def winSd_lp_loader(file_path):
    # Dictionary to store the extracted values
    extracted_data = {}

    # Define the keys you want to extract
    keys = {
        "Le": "Le",
        "Re": "Re",
        "BL": "Bl",  # Case-insensitive key matching for Bl
        "Sd": "Sd",
        "Cms": "Cms",
        "Mms": "Mms",
        "Rms": "Rms"
    }

    # Open the file and read line by line
    with open(file_path, 'r') as file:
        for line in file:
            # Split each line by '=' and ensure valid format
            if '=' in line:
                key, value = line.strip().split('=')
                key = key.strip()  # Clean up the key
                value = value.strip()  # Clean up the value

                # Check if the key is one of the desired ones and not already extracted
                if key in keys and keys[key] not in extracted_data:
                    extracted_data[keys[key]] = float(value)

    return extracted_data


def klippel_lp_loader(file_path):
    # Dictionary to store the extracted values
    extracted_data = {}

    # Define the keys you want to extract
    keys = {
        "Re": "Re",
        "Le": "Le",
        "Bl": "Bl",
        "Sd": "Sd",
        "Cms": "Cms",
        "Mms": "Mms",
        "Rms": "Rms"
    }

    # Open the file and read line by line
    with open(file_path, 'r') as file:
        for line in file:
            # Split the line by tab or spaces if it contains a parameter
            if any(key in line for key in keys):
                parts = line.strip().split()  # Split the line by whitespace

                if len(parts) > 1:
                    key = parts[0]  # First part of the line is the key
                    value = parts[1]  # Second part of the line is the value

                    # Check if the key is one of the desired ones and not already extracted
                    if key in keys and keys[key] not in extracted_data:
                        extracted_data[keys[key]] = float(value)

    return extracted_data





