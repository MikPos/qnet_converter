rule [
    left [
        edge [ source 2 target 1 label "=" ]
    ]
    context [
        node [ id 1 label "C" ]
        node [ id 2 label "C" ]
        node [ id 3 label "O" ]
    ]
    right [
        edge [ source 2 target 1 label "-" ]
        edge [ source 2 target 3 label "-" ]
    ]
]
