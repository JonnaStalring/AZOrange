import os


#APWS = "192.168.100.238:8082"

# Define the ADMET Predictor endpoints for individual calculation
APENDPOINTS = ["logP", "logD", "MDCK", "Peff", "pKa", "RuleOf3", "RuleOf5", "Sp", "Acidic_pKa", "Acidic_pKa_74prox", "Basic_pKa", "Basic_pKa_74prox", "Mixed_pKa", "Mixed_pKa_74prox"]
# The AP endpoints included in the AllAPendpoints endpoint
ALLAPENDPOINTSLIST = ["logP", "logD", "MDCK", "Peff", "RuleOf3", "RuleOf5", "Sp", "Acidic_pKa", "Acidic_pKa_74prox", "Basic_pKa", "Basic_pKa_74prox", "Mixed_pKa", "Mixed_pKa_74prox", "pKa_mostAcidic", "pKa_mostBasic"]

# The maximum number of molecules that can be processed for an AP endpoints
NAPMOL = 10 
NAPMOLBATCH = 1500

# Define the directories
CHEMICSROOTPATH = os.environ['CHEMICSROOTPATH']
CHEMICSMODELCODEPATH = os.path.join(CHEMICSROOTPATH, "ChemicsEndpoints/Endpoints")
CHEMICSMODELDIR = os.path.join(CHEMICSROOTPATH, "ChemicsModelDir")

#CHEMICSTIMEOUT = 300
CHEMICSTIMEOUT = 5  # Keep while feedback while waiting for license is not implemented

ERRORCODE = "Error"
FINISHEDCODE = "Finished"
S_ERRORCODE = "NotStarted"
APLICENSEBUSY = "ADMET Predictor License busy"

# Retrun codes from getStatus
QUEUEDCODE = "Queued"
RUNNINGCODE = "Running"
SUCCESSCODE = "Completed: All molecules predicted successfully"
PARTIALSUCCESSCODE = "Incomplete: Some molecules could not be predicted. Please see the 'Calculation status' column. In case of errors, please report to Helpdesk providing the information in this box (Copy Summary To Clipboard)."
FAILEDCODE = "TASK FAILED: In case of errors, please report to Helpdesk providing the information in this box (Copy Summary To Clipboard)."


D360ENDPOINTS = [
{"HeavyAtomCount": {"unit": "", "version": "RDK12.12.1"}},
{"AtomCount": {"unit": "", "version": "RDK12.12.1"}},
{"MolWt": {"unit": "g/mol", "version": "RDK12.12.1"}},
{"BondCount": {"unit": "", "version": "RDK12.12.1"}},
{"FluorineCount": {"unit": "", "version": "RDK12.12.1"}},
{"HalogenCount": {"unit": "", "version": "RDK12.12.1"}},
{"CarbonCount": {"unit": "", "version": "RDK12.12.1"}},
{"PhosphorusCount": {"unit": "", "version": "RDK12.12.1"}},
{"ChlorineCount": {"unit": "", "version": "RDK12.12.1"}},
{"SulfurCount": {"unit": "", "version": "RDK12.12.1"}},
{"NitrogenCount": {"unit": "", "version": "RDK12.12.1"}},
{"OxygenCount": {"unit": "", "version": "RDK12.12.1"}},
{"SMILES": {"unit": "", "version": "RDK12.12.1"}},
{"RotatableBondsCount": {"unit": "", "version": "RDK12.12.1"}},
{"RingCount": {"unit": "", "version": "RDK12.12.1"}},
{"HDonorsCount": {"unit": "", "version": "RDK12.12.1"}},
{"TPSA": {"unit": "Angstrom^2", "version": "RDK12.12.1"}},
{"HAcceptorsCount": {"unit": "", "version": "RDK12.12.1"}},
{"logP": {"unit": "", "version": "AP7.1"}},
{"logD": {"unit": "", "version": "AP7.1"}},
{"MDCK": {"unit": "cm/s*10^7", "version": "AP7.1"}},
{"Peff": {"unit": "cm/s*10^4", "version": "AP7.1"}},
{"RuleOf3": {"unit": "", "version": "AP7.1"}},
{"RuleOf5": {"unit": "", "version": "AP7.1"}},
{"Sp": {"unit": "mg/mL", "version": "AP7.1"}},
{"Acidic_pKa": {"unit": "", "version": "AP7.1"}},
{"Acidic_pKa_74prox": {"unit": "", "version": "AP7.1"}},
{"Basic_pKa": {"unit": "", "version": "AP7.1"}},
{"Basic_pKa_74prox": {"unit": "", "version": "AP7.1"}},
{"Mixed_pKa": {"unit": "", "version": "AP7.1"}},
{"Mixed_pKa_74prox": {"unit": "", "version": "AP7.1"}},
{"pKa_mostBasic": {"unit": "", "version": "AP7.1"}},
{"pKa_mostAcidic": {"unit": "", "version": "AP7.1"}}
]


# Define the endpoints and their hierachy for D360
D360ENDPOINTHIERARCHY = {
  "type": "folder",
  "name": "Chemics Service",
  "displayname": "Chemics Service",
  "version": None,
  "outputProperties": None,
  "subFolders":   [
        # General
        {
        "type": "folder",
        "name": "General",
        "displayname": "General",
        "version": None,
        "outputProperties": None,
        "subFolders":       [
            {
                "type": "task",
                "name": "SMILES",
                "displayname": "SMILES_RDK12.12.1",
                "version": "RDK12.12.1",
                "outputProperties": [ 
                    {
                    "descriptorId": "prediction",
                    "displayname": "SMILES_RDK12.12.1",
                    "dataType": "string",
                    "unit": None,
                    "selectedByDefault": True
                    },
                    {
                    "descriptorId": "calcStatus",
                    "displayname": "Calculation status SMILES",
                    "dataType": "string",
                    "unit": None,
                    "selectedByDefault": True
                    }
                    ],  # End outputProperties
                    "subFolders": None 
            },    # Finish task
            {
                "type": "task",
                "name": "MolWt",
                "displayname": "MolWt_RDK12.12.1",
                "version": "RDK12.12.1",
                "outputProperties": [
                    {
                    "descriptorId": "prediction",
                    "displayname": "MolWt_RDK12.12.1 (g/mol)",
                    "dataType": "float",
                    "unit": "g/mol",
                    "selectedByDefault": True
                    },
                    {
                    "descriptorId": "calcStatus",
                    "displayname": "Calculation status MolWt",
                    "dataType": "string",
                    "unit": None,
                    "selectedByDefault": True
                    }
                    ],  # End outputProperties
                    "subFolders": None
            }    # Finish task
            ]    # Finish General subfolders
        },
        # Counts
        {
        "type": "folder",
        "name": "Counts",
        "displayname": "Counts",
        "version": None,
        "outputProperties": None,
        "subFolders":       [
            {
                "type": "task",
                "name": "HeavyAtomCount",
                "displayname": "HeavyAtomCount_RDK12.12.1",
                "version": "RDK12.12.1",
                "outputProperties": [
                    {
                    "descriptorId": "prediction",
                    "displayname": "HeavyAtomCount_RDK12.12.1",
                    "dataType": "int",
                    "unit": None,
                    "selectedByDefault": True
                    },
                    {
                    "descriptorId": "calcStatus",
                    "displayname": "Calculation status HeavyAtomCount",
                    "dataType": "string",
                    "unit": None,
                    "selectedByDefault": True
                    }
                    ],  # End outputProperties
                    "subFolders": None
            },    # Finish task
            {
                "type": "task",
                "name": "BondCount",
                "displayname": "BondCount_RDK12.12.1",
                "version": "RDK12.12.1",
                "outputProperties": [
                    {
                    "descriptorId": "prediction",
                    "displayname": "BondCount_RDK12.12.1",
                    "dataType": "int",
                    "unit": None,
                    "selectedByDefault": True
                    },
                    {
                    "descriptorId": "calcStatus",
                    "displayname": "Calculation status BondCount",
                    "dataType": "string",
                    "unit": None,
                    "selectedByDefault": True
                    }
                    ],  # End outputProperties
                    "subFolders": None
            },    # Finish task
            {
                "type": "task",
                "name": "FluorineCount",
                "displayname": "FluorineCount_RDK12.12.1",
                "version": "RDK12.12.1",
                "outputProperties": [
                    {
                    "descriptorId": "prediction",
                    "displayname": "FluorineCount_RDK12.12.1",
                    "dataType": "int",
                    "unit": None,
                    "selectedByDefault": True
                    },
                    {
                    "descriptorId": "calcStatus",
                    "displayname": "Calculation status FluorineCount",
                    "dataType": "string",
                    "unit": None,
                    "selectedByDefault": True
                    }
                    ],  # End outputProperties
                    "subFolders": None
            },    # Finish task
            {
                "type": "task",
                "name": "HalogenCount",
                "displayname": "HalogenCount_RDK12.12.1",
                "version": "RDK12.12.1",
                "outputProperties": [
                    {
                    "descriptorId": "prediction",
                    "displayname": "HalogenCount_RDK12.12.1",
                    "dataType": "int",
                    "unit": None,
                    "selectedByDefault": True
                    },
                    {
                    "descriptorId": "calcStatus",
                    "displayname": "Calculation status HalogenCount",
                    "dataType": "string",
                    "unit": None,
                    "selectedByDefault": True
                    }
                    ],  # End outputProperties
                    "subFolders": None
            },    # Finish task
            {
                "type": "task",
                "name": "CarbonCount",
                "displayname": "CarbonCount_RDK12.12.1",
                "version": "RDK12.12.1",
                "outputProperties": [
                    {
                    "descriptorId": "prediction",
                    "displayname": "CarbonCount_RDK12.12.1",
                    "dataType": "int",
                    "unit": None,
                    "selectedByDefault": True
                    },
                    {
                    "descriptorId": "calcStatus",
                    "displayname": "Calculation status CarbonCount",
                    "dataType": "string",
                    "unit": None,
                    "selectedByDefault": True
                    }
                    ],  # End outputProperties
                    "subFolders": None
            },    # Finish task
            {
                "type": "task",
                "name": "PhosphorusCount",
                "displayname": "PhosphorusCount_RDK12.12.1",
                "version": "RDK12.12.1",
                "outputProperties": [
                    {
                    "descriptorId": "prediction",
                    "displayname": "PhosphorusCount_RDK12.12.1",
                    "dataType": "int",
                    "unit": None,
                    "selectedByDefault": True
                    },
                    {
                    "descriptorId": "calcStatus",
                    "displayname": "Calculation status PhosphorusCount",
                    "dataType": "string",
                    "unit": None,
                    "selectedByDefault": True
                    }
                    ],  # End outputProperties
                    "subFolders": None
            },    # Finish task
            {
                "type": "task",
                "name": "ChlorineCount",
                "displayname": "ChlorineCount_RDK12.12.1",
                "version": "RDK12.12.1",
                "outputProperties": [
                    {
                    "descriptorId": "prediction",
                    "displayname": "ChlorineCount_RDK12.12.1",
                    "dataType": "int",
                    "unit": None,
                    "selectedByDefault": True
                    },
                    {
                    "descriptorId": "calcStatus",
                    "displayname": "Calculation status ChlorineCount",
                    "dataType": "string",
                    "unit": None,
                    "selectedByDefault": True
                    }
                    ],  # End outputProperties
                    "subFolders": None
            },    # Finish task
            {
                "type": "task",
                "name": "OxygenCount",
                "displayname": "OxygenCount_RDK12.12.1",
                "version": "RDK12.12.1",
                "outputProperties": [
                    {
                    "descriptorId": "prediction",
                    "displayname": "OxygenCount_RDK12.12.1",
                    "dataType": "int",
                    "unit": None,
                    "selectedByDefault": True
                    },
                    {
                    "descriptorId": "calcStatus",
                    "displayname": "Calculation status OxygenCount",
                    "dataType": "string",
                    "unit": None,
                    "selectedByDefault": True
                    }
                    ],  # End outputProperties
                    "subFolders": None
            },    # Finish task
            {
                "type": "task",
                "name": "TPSA",
                "displayname": "TPSA_RDK12.12.1",
                "version": "RDK12.12.1",
                "outputProperties": [
                    {
                    "descriptorId": "prediction",
                    "displayname": "TPSA_RDK12.12.1",
                    "dataType": "int",
                    "unit": None,
                    "selectedByDefault": True
                    },
                    {
                    "descriptorId": "calcStatus",
                    "displayname": "Calculation status TPSA",
                    "dataType": "string",
                    "unit": None,
                    "selectedByDefault": True
                    }
                    ],  # End outputProperties
                    "subFolders": None
            },    # Finish task
            {
                "type": "task",
                "name": "HAcceptorsCount",
                "displayname": "HAcceptorsCount_RDK12.12.1",
                "version": "RDK12.12.1",
                "outputProperties": [
                    {
                    "descriptorId": "prediction",
                    "displayname": "HAcceptorsCount_RDK12.12.1",
                    "dataType": "int",
                    "unit": None,
                    "selectedByDefault": True
                    },
                    {
                    "descriptorId": "calcStatus",
                    "displayname": "Calculation status HAcceptorsCount",
                    "dataType": "string",
                    "unit": None,
                    "selectedByDefault": True
                    }
                    ],  # End outputProperties
                    "subFolders": None
            },    # Finish task
            {
                "type": "task",
                "name": "HDonorsCount",
                "displayname": "HDonorsCount_RDK12.12.1",
                "version": "RDK12.12.1",
                "outputProperties": [
                    {
                    "descriptorId": "prediction",
                    "displayname": "HDonorsCount_RDK12.12.1",
                    "dataType": "int",
                    "unit": None,
                    "selectedByDefault": True
                    },
                    {
                    "descriptorId": "calcStatus",
                    "displayname": "Calculation status HDonorsCount",
                    "dataType": "string",
                    "unit": None,
                    "selectedByDefault": True
                    }
                    ],  # End outputProperties
                    "subFolders": None
            },    # Finish task
            {
                "type": "task",
                "name": "RingCount",
                "displayname": "RingCount_RDK12.12.1",
                "version": "RDK12.12.1",
                "outputProperties": [
                    {
                    "descriptorId": "prediction",
                    "displayname": "RingCount_RDK12.12.1",
                    "dataType": "int",
                    "unit": None,
                    "selectedByDefault": True
                    },
                    {
                    "descriptorId": "calcStatus",
                    "displayname": "Calculation status RingCount",
                    "dataType": "string",
                    "unit": None,
                    "selectedByDefault": True
                    }
                    ],  # End outputProperties
                    "subFolders": None
            },    # Finish task
            {
                "type": "task",
                "name": "RotatableBondsCount",
                "displayname": "RotatableBondsCount_RDK12.12.1",
                "version": "RDK12.12.1",
                "outputProperties": [
                    {
                    "descriptorId": "prediction",
                    "displayname": "RotatableBondsCount_RDK12.12.1",
                    "dataType": "int",
                    "unit": None,
                    "selectedByDefault": True
                    },
                    {
                    "descriptorId": "calcStatus",
                    "displayname": "Calculation status RotatableBondsCount",
                    "dataType": "string",
                    "unit": None,
                    "selectedByDefault": True
                    }
                    ],  # End outputProperties
                    "subFolders": None
            },    # Finish task
            {
                "type": "task",
                "name": "NitrogenCount",
                "displayname": "NitrogenCount_RDK12.12.1",
                "version": "RDK12.12.1",
                "outputProperties": [
                    {
                    "descriptorId": "prediction",
                    "displayname": "NitrogenCount_RDK12.12.1",
                    "dataType": "int",
                    "unit": None,
                    "selectedByDefault": True
                    },
                    {
                    "descriptorId": "calcStatus",
                    "displayname": "Calculation status NitrogenCount",
                    "dataType": "string",
                    "unit": None,
                    "selectedByDefault": True
                    }
                    ],  # End outputProperties
                    "subFolders": None
            },    # Finish task
            {
                "type": "task",
                "name": "SulfurCount",
                "displayname": "SulfurCount_RDK12.12.1",
                "version": "RDK12.12.1",
                "outputProperties": [
                    {
                    "descriptorId": "prediction",
                    "displayname": "SulfurCount_RDK12.12.1",
                    "dataType": "int",
                    "unit": None,
                    "selectedByDefault": True
                    },
                    {
                    "descriptorId": "calcStatus",
                    "displayname": "Calculation status SulfurCount",
                    "dataType": "string",
                    "unit": None,
                    "selectedByDefault": True
                    }
                    ],  # End outputProperties
                    "subFolders": None
            },    # Finish task
            {
                "type": "task",
                "name": "AtomCount",
                "displayname": "AtomCount_RDK12.12.1",
                "version": "RDK12.12.1",
                "outputProperties": [
                    {
                    "descriptorId": "prediction",
                    "displayname": "AtomCount_RDK12.12.1",
                    "dataType": "int",
                    "unit": None,
                    "selectedByDefault": True
                    },
                    {
                    "descriptorId": "calcStatus",
                    "displayname": "Calculation status AtomCount",
                    "dataType": "string",
                    "unit": None,
                    "selectedByDefault": True
                    }
                    ],  # End outputProperties
                    "subFolders": None
            }    # Last Counts object
            ]    # Finish Counts subfolders
        }, # Finish Counts
        # AllAPEndpoints
        {
        "type": "folder",
        "name": "AllAPendpoints",
        "displayname": "AllAPendpoints",
        "version": None,
        "outputProperties": None,
        "subFolders":    [
            {
                "type": "task",
                "name": "AllAPendpoints",
                "displayname": "AllAPendpoints_AP7.1",
                "version": "AP7.1",
                "outputProperties": [
                    {
                    "descriptorId": "prediction",
                    "displayname": "AllAPendpoints_AP7.1",
                    "dataType": "string",
                    "unit": None,
                    "selectedByDefault": True
                    },
                    {
                    "descriptorId": "calcStatus",
                    "displayname": "Calculation status AllAPendpoints",
                    "dataType": "string",
                    "unit": None,
                    "selectedByDefault": True
                    }
                    ],  # End outputProperties
                    "subFolders": None
            },    # Finish task
            ]    # Finish AllAPndpoints subfolders
        },
        # PhysChem
        {
        "type": "folder",
        "name": "PhysChem",
        "displayname": "PhysChem",
        "version": None,
        "outputProperties": None,
        "subFolders":       [
            {
            "type": "folder",
            "name": "Lipophilicity",
            "displayname": "Lipophilicity",
            "version": None,
            "outputProperties": None,
            "subFolders":   None,
            },   # Lipophilicity folder
            {
            "type": "folder",
            "name": "Solubility",
            "displayname": "Solubility",
            "version": None,
            "outputProperties": None,
            "subFolders":    None,
             }, #Solubility folder
             {
             "type": "folder",
             "name": "pKa",
             "displayname": "pKa",
             "version": None,
             "outputProperties": None,
             "subFolders":       None,
             }, # pKa folder
             {
             "type": "folder",
             "name": "Rules",
             "displayname": "Rules",
             "version": None,
             "outputProperties": None,
             "subFolders":   None
             }   # Rules
            ] # PhysChem sub folders
        },  # PhysChem folder
        {
        "type": "folder",
        "name": "DMPK",
        "displayname": "DMPK",
        "version": None,
        "outputProperties": None,
        "subFolders":       [
                    {
                    "type": "folder",
                    "name": "Permeability",
                    "displayname": "Permeability",
                    "version": None,
                    "outputProperties": None,
                    "subFolders":   None,
                    } # Permeability folder
            ]  # DMPK subfolders
        }  # DMPK folder
        ]  # Finish Chemics Service subfolders
}

#"General":
#    [
#    {"SMILES": {"unit": "", "version": "RDK12.12.1"}}
#    ],
#"Counts":
#    [
#    {"HeavyAtomCount": {"unit": "", "version": "RDK12.12.1"}},
#    {"AtomCount": {"unit": "", "version": "RDK12.12.1"}},
#    {"BondCount": {"unit": "", "version": "RDK12.12.1"}},
#    {"FluorineCount": {"unit": "", "version": "RDK12.12.1"}},
#    {"HalogenCount": {"unit": "", "version": "RDK12.12.1"}},
#    {"CarbonCount": {"unit": "", "version": "RDK12.12.1"}},
#    {"PhosphorusCount": {"unit": "", "version": "RDK12.12.1"}},
#    {"ChlorineCount": {"unit": "", "version": "RDK12.12.1"}},
#    {"SulfurCount": {"unit": "", "version": "RDK12.12.1"}},
#    {"NitrogenCount": {"unit": "", "version": "RDK12.12.1"}},
#    {"OxygenCount": {"unit": "", "version": "RDK12.12.1"}},
#    {"RotatableBondsCount": {"unit": "", "version": "RDK12.12.1"}},
#    {"RingCount": {"unit": "", "version": "RDK12.12.1"}},
#    {"HDonorsCount": {"unit": "", "version": "RDK12.12.1"}},
#    {"HAcceptorsCount": {"unit": "", "version": "RDK12.12.1"}},
#    {"TPSA": {"unit": "Angstrom^2", "version": "RDK12.12.1"}}
#    ],
#"PhysChem":
#    {
#    "Lipophilicity":
#        [
#        {"logP": {"unit": "", "version": "AP7.1"}},
#        {"logD": {"unit": "", "version": "AP7.1"}}
#        ],
#    "Solubility":
#        [
#        {"Sp": {"unit": "mg/mL", "version": "AP7.1"}}
#        ],
#    "pKa":
#        [
#        {"Acidic_pKa": {"unit": "", "version": "AP7.1"}},
#        {"Acidic_pKa_74prox": {"unit": "", "version": "AP7.1"}},
#        {"Basic_pKa": {"unit": "", "version": "AP7.1"}},
#        {"Basic_pKa_74prox": {"unit": "", "version": "AP7.1"}},
#        {"Mixed_pKa": {"unit": "", "version": "AP7.1"}},
#        {"Mixed_pKa_74prox": {"unit": "", "version": "AP7.1"}}
#        ],
#    },
#"Rules":
#    [
#    {"RuleOf3": {"unit": "", "version": "AP7.1"}},
#    {"RuleOf5": {"unit": "", "version": "AP7.1"}}
#    ],
#"DMPK":
#    {
#    "Permeability":
#        [
#        {"MDCK": {"unit": "cm/s*10^7", "version": "AP7.1"}},
#        {"Peff": {"unit": "cm/s*10^4", "version": "AP7.1"}}
#        ],
#    "Metabolism":
#        [
#        ],
#    },
#"Safety":
#    [
#    ]
#}

