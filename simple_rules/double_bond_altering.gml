rule [
    left [
        node [ id 5 label "O." ]
        node [ id 1 label "C" ]
        edge [ source 1 target 2 label "=" ]
        edge [ source 2 target 3 label "-" ]
        edge [ source 3 target 4 label "-" ]
    ]
    context [
        node [ id 2 label "C" ]
        node [ id 3 label "O" ]
        node [ id 4 label "H" ]
    ]
    right [
        node [ id 5 label "O" ]
        node [ id 1 label "C." ]
        edge [ source 5 target 4 label "-" ]
        edge [ source 1 target 2 label "-" ]
        edge [ source 2 target 3 label "=" ]
    ]
]
