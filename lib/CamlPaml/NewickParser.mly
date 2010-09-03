%token <string> LABEL
%token <float> BRANCHLEN
%token COLON LPAREN RPAREN COMMA SEMICOLON EOF
%start parse
%type <Newick.t> parse
%%
parse:
	| node SEMICOLON								{ $1 }
	| node EOF										{ $1 }
;
node:
	| label											{ Newick.Node ([],fst $1,snd $1) }
	| LPAREN subtrees RPAREN label					{ Newick.Node ($2,fst $4,snd $4) }
	| LPAREN subtrees RPAREN						{ Newick.Node ($2,"",None) }
;
label:
	| LABEL											{ ($1,None) }
	| LABEL COLON BRANCHLEN							{ ($1,Some $3) }
	| COLON BRANCHLEN								{ ("",Some $2) }
;
subtrees:
	| node COMMA subtrees							{ $1 :: $3 }
	| node											{ [$1] }
;
%%
