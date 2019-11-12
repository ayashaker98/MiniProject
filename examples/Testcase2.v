module test2();
input [3:0]s;
input [2:0]L;
output n;
AND2X1 AND2X1_1 ( .A(s[0]), .B(s[1]), .Y(c1) );
OAI21X1 OAI21X1_1 ( .A(s[1]), .B(s[2]), .C(c1), .Y(c2) );
AND2X1 AND2X1_2 ( .A(s[2]), .B(c1), .Y(c3) );
INVX1 INVX1_1 ( .A(c3), .Y(c4) );
AND2X1 AND2X1_3 ( .A(s[0]), .B(s[1]), .Y(c5) );
AND2X1 AND2X1_4 ( .A(s[0]), .B(s[1]), .Y(c6) );
AND2X1 AND2X1_5 ( .A(s[1]), .B(s[2]), .Y(c7) );
NOR2X1  NOR2X1_1 ( .A(c7), .B(c4), .Y(c8) );
AND2X1 AND2X1_6 ( .A(L[0]), .B(L[1]), .Y(c9) );
AND2X1 AND2X1_7 ( .A(L[2]), .B(c9), .Y(c10) );
AND2X1 AND2X1_8 ( .A(c10), .B(c9), .Y(c11) );
OR2X1  OR2X1_1 ( .A(c11), .B(L[2]), .Y(n) );
endmodule