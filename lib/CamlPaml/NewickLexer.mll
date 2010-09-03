{
open NewickParser
}

rule token = parse
	| [' ' '\t' '\r' '\n' ';']							{ token lexbuf }
	| '('												{ LPAREN }
	| ')'												{ RPAREN }
	| ','												{ COMMA }
	| ':'												{ COLON }
	| ['0'-'9' '.']+ as lxm								{ BRANCHLEN(float_of_string lxm) }
	| ['A'-'Z' 'a'-'z' '0'-'9' '_' '.']+ as lxm			{ LABEL(lxm) }
	| ';'												{ SEMICOLON }
	| eof												{ EOF }
