# Comparação de rotinas para resolução de equações do tipo TOV

Comparamos as diferentes linguagens na solução de problemas do tipo TOV (solução das equações de Tolman-Oppenheimer-Volkoff), mantendo o mesmo algorítimo para todos.

## Dependências

O código em Julia necessita de algumas bibliotecas instaladas a partir do próprio diretório `julia` e REPL com:

```
$ julia --project=.
julia> ]
pkg> add https://github.com/SunkPlane29/tov.git
pkg> instantiate
```

---

O código em Fortran necessita das bibliotecas Blas e Lapack que, caso não presentes, podem ser instaladas no Ubuntu (e similares) com:

```
$ sudo apt-get install libblas-dev liblapack-dev
```

---

O código em Python necessita das bibliotecas:
- Numpy
- Matplotlib
- Pandas

Elas podem ser instaladas com:

```
$ pip install numpy matplotlib pandas
```

---

Ao rodar o código Júlia pode haver algum problema com a dependência TOV. Basta instalar ela novamente no ambiente:

```
$ julia --project=.
julia> ] add https://github.com/SunkPlane29/tov.git
```
