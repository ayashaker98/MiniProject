module test1(s,r,T,f):
    input [3:0]s;
    input [3:0]r;
    output T;
    output f;

BUFX2 BUFX2_1 ( .A(s[0]), .Y(c1) );
BUFX2 BUFX2_2 ( .A(s[1]), .Y(c2) );
BUFX2 BUFX2_3 ( .A(s[2]), .Y(c3) );
BUFX2 BUFX2_4 ( .A(s[3]), .Y(c4) );
AND2X1 AND2X1_1 ( .A(c3), .B(s[3]), .Y(L1) );
NAND2X1 NAND2X1_1 ( .A(c1), .B(s[3]), .Y(L2) );
NAND2X1 NAND2X1_2 ( .A(c1), .B(s[2]), .Y(L3) );
AND2X1 AND2X1_2 ( .A(c1), .B(s[3]), .Y(L4) );
INVX1 INVX1_1 ( .A(L4), .Y(Y5) );
INVX1 INVX1_2 ( .A(L4), .Y(L6) );
OAI21X1 OAI21X1_3 ( .A(c2), .B(r[3]), .C(r[2]), .Y(T) );
OR2X1 OR2X1_1 ( .A(L3), .B(L1), .Y(n));
OR2X1 OR2X1_2 ( .A(n), .B(r[0]), .Y(T));
OR2X1 OR2X1_3 ( .A(L4), .B(c4), .Y(f));

endmodule