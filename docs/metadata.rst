Metadata format
===============

Project metadata format::

    {
        "HomologSet" : [
            {
                "id" : "XX00000000",
                "sequence" : "AGAMAMGATKLLSMA",
                "orgname" : "human",
            },
            {
                "id" : "XX00000001",
                "sequence" : "AGAKKLGATKLLSMA",
                "orgname" : "human",
            }
        ],
        "Alignment" : [
            {
                "XX00000000" : "AG---KKLGATKL"
            },
            {
                "XX00000001" : ""
            },
        ],
        ""
    }

HomologSet metadata::

    {"HomologSet" : [
            {
                "id" : "XX00000000",
                "sequence" : "AGAMAMGATKLLSMA",
                "orgname" : "human",
            },
            {
                "id" : "XX00000001",
                "sequence" : "AGAKKLGATKLLSMA",
                "orgname" : "human",
            }
        ]
    }

Alignment metadata format::
    {"Alignment" : }
