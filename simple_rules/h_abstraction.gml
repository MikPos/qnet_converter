rule [
    left [
        node [ id 3 label "O." ]
        node [ id 1 label "O" ]
        edge [ source 2 target 1 label "-" ]
    ]
    context [
        node [ id 2 label "H" ]
    ]
    right [
        node [ id 3 label "O" ]
        node [ id 1 label "O." ]
        edge [ source 2 target 3 label "-" ]
    ]
]
