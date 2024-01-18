#include <stdio.h>
#include <stdlib.h>
#include <GL/glut.h>
#include <cmath>
#include <fstream>
#include <iostream>
#define MAX_PONTOS 1000
#define MAX_RETAS 1000
#define MAX_POLYGON_POINTS 100
#define MAX_VERTICES 100

using namespace std;

typedef struct { //Definição da estrutura Ponto
    float x;
    float y;
    float cor[3];
} Ponto;

typedef struct { //Definição da estrutura Reta
    Ponto inicio;
    Ponto fim;
    float cor[3];
} Reta;

typedef struct { //Definição da estrutura Poligono
    int quantidade_vertices;
    Ponto vertices[MAX_VERTICES];
    float cor[3];
} Poligono;

// Variáveis globais
int modo_desenho = 0;
int quantidade_pontos = 0;
Ponto pontos[MAX_PONTOS];
int quantidade_retas = 0;
Reta retas[MAX_RETAS];
int quantidade_poligonos = 0;
Poligono poligonos[MAX_PONTOS];
Ponto pontosPoligono[MAX_POLYGON_POINTS];
int quantidadePontosPoligono = 0;
int desenhaPoligono = 0;
int transformacao = 0;
Ponto pontoInicialTransformacao;
int objetoatual = 0;

// Funções de manipulação de estruturas
Ponto criarPonto(float x, float y, float cor[3]) {
    Ponto p;
    p.x = x;
    p.y = y;
    for (int i = 0; i < 3; i++) {
        p.cor[i] = cor[i];
    }
    return p;
}

Reta criarReta(float x1, float y1, float x2, float y2, float cor[3]) {
    Reta r;
    r.inicio = criarPonto(x1, y1, cor);
    r.fim = criarPonto(x2, y2, cor);
    for (int i = 0; i < 3; i++) {
        r.cor[i] = cor[i];
    }
    return r;
}

Poligono criarPoligono(Ponto pontos[], int numPontos, float cor[3]) {
    Poligono p;
    p.quantidade_vertices = numPontos;
    for (int i = 0; i < numPontos; i++) {
        p.vertices[i] = pontos[i];
    }
    for (int i = 0; i < 3; i++) {
        p.cor[i] = cor[i];
    }
    return p;
} //Fim das Funções de manipulação de estruturas

// Função para multiplicação de matrizes
void multiplicarMatrizes(float result[3][1], float matriz1[3][3], float matriz2[3][1]) {
    for (int i = 0; i < 3; i++) {
        result[i][0] = 0;
        for (int k = 0; k < 3; k++) {
            result[i][0] += matriz1[i][k] * matriz2[k][0];
        }
    }
}

// Função para copiar matrizes
void copiarMatrizes(float matriz1[3][3], float matriz2[3][3]) {
    for (int i = 0; i < 3; i++) {
        for (int k = 0; k < 3; k++) {
            matriz1[i][k] = matriz2[i][k];
        }
    }
}

// Funções para adicionar elementos
void adicionarPonto(float x, float y, float cor[3]) {
    pontos[quantidade_pontos++] = criarPonto(x, y, cor);
}

void adicionarReta(float x1, float y1, float x2, float y2, float cor[3]) {
    retas[quantidade_retas++] = criarReta(x1, y1, x2, y2, cor);
}

void adicionarPoligono(Ponto pontos[], int numPontos, float cor[3]) {
    poligonos[quantidade_poligonos++] = criarPoligono(pontos, numPontos, cor);
}// Fim das Funções para adicionar elementos


// Funções para desenhar elementos
void desenharPontos() {
    glPointSize(5.0);
    glBegin(GL_POINTS);
    for (int i = 0; i < quantidade_pontos; i++) {
        glColor3fv(pontos[i].cor);
        glVertex2f(pontos[i].x, pontos[i].y);
    }
    glEnd();
}

void desenharRetas() {
    glLineWidth(2.0);
    glBegin(GL_LINES);
    for (int i = 0; i < quantidade_retas; i++) {
        glColor3fv(retas[i].inicio.cor);
        glVertex2f(retas[i].inicio.x, retas[i].inicio.y);
        glColor3fv(retas[i].fim.cor);
        glVertex2f(retas[i].fim.x, retas[i].fim.y);
    }
    glEnd();
}

void desenharPoligonos() {
    for (int i = 0; i < quantidade_poligonos; i++) {
        glBegin(GL_POLYGON);
        glColor3fv(poligonos[i].cor);
        for (int j = 0; j < poligonos[i].quantidade_vertices; j++) {
            glVertex2f(poligonos[i].vertices[j].x, poligonos[i].vertices[j].y);
        }
        glEnd();
    }
}// Fim das Funções para desenhar elementos
// Função que verifica se um ponto (mx, my) está dentro da margem de um ponto alvo (px, py) com tolerância t.
int selecPonto(float px, float py, float mx, float my, int t) {
    float margem = 0.01; //Precisei usar essa margem pois por algum motivo nao consegui fazer funcionar precisamente só com t.
    if (mx <= px + t + margem && mx >= px - t - margem) {
        if (my <= py + t + margem && my >= py - t - margem) {
            return 1;// O ponto (mx, my) está dentro da margem do ponto alvo (px, py).
        }
    }
    return 0;// O ponto (mx, my) está fora da margem do ponto alvo (px, py).
}
// Função que percorre um conjunto de pontos e retorna o índice do primeiro ponto dentro da tolerância t.
int retornaPonto(float mx, float my, int tolerancia) {
    for (int i = 0; i < quantidade_pontos; i++) {
        if (selecPonto(pontos[i].x, pontos[i].y, mx, my, tolerancia)) {
            return i; // Retorna o índice do primeiro ponto dentro da tolerância.
        }
    }
    return -1; // Nenhum ponto encontrado dentro da tolerância.
}
// Função que verifica um ponto dentro da tolerância em relação às coordenadas (mx, my) e o remove se encontrado.
int pontoSelecionado(float mx, float my, int tolerancia) {
    for (int i = 0; i < quantidade_pontos; i++) {
        if (selecPonto(pontos[i].x, pontos[i].y, mx, my, tolerancia)) { // Remove o ponto encontrado e ajusta vetor.
            for (int j = i; j < quantidade_pontos - 1; j++) {
                pontos[j] = pontos[j + 1];
            }
            quantidade_pontos--; // Atualiza a quantidade de pontos após a remoção.
            return 1; // Retorna 1 indicando que um ponto foi removido.
        }
    }
    return 0;// Retorna 0 indicando que nenhum ponto foi removido.
}

const int NUM_BITS = 4;
// Função que codifica um ponto com base nas comparações com os limites de coordenadas: xmin, xmax, ymin, ymax.
int codificarPonto(Ponto ponto, float xmin, float xmax, float ymin, float ymax) {
    int codigo = 0;
    // Fiz seguindo as anotações dos slides vistos em aula.
    if (ponto.x < xmin) codigo |= 1; // Bit 0 representa a posição à esquerda do limite xmin.
    if (ponto.x > xmax) codigo |= 2; // Bit 1 representa a posição à direita do limite xmax.
    if (ponto.y < ymin) codigo |= 4; // Bit 2 representa a posição abaixo do limite ymin.
    if (ponto.y > ymax) codigo |= 8; // Bit 3 representa a posição acima do limite ymax.

    return codigo;
}
// Função que verifica se há pelo menos um bit comum entre dois códigos.
int verificaBit(int codigo1, int codigo2) {
    return (codigo1 & codigo2) != 0;
}
// Função de seleção da reta.
int selecionarReta(Reta reta, float tolerancia, float mx, float my) {
    //Limites da janela com base na tolerância.
    float xmin = mx - tolerancia;
    float xmax = mx + tolerancia;
    float ymin = my - tolerancia;
    float ymax = my + tolerancia;
    // Codifica os pontos iniciais da reta em relação à janela.
    int codigoP0 = codificarPonto(reta.inicio, xmin, xmax, ymin, ymax);//ponto inicial.
    int codigoP1 = codificarPonto(reta.fim, xmin, xmax, ymin, ymax);//ponto final.
    // Verifica se há interseção inicial.
    if (verificaBit(codigoP0, codigoP1))
        return 0;
    // Verifica se ambos os pontos estão dentro da janela inicialmente.
    if (codigoP0 == 0 && codigoP1 == 0)
        return 1;

    while (1) {// Aplicando a logica do algoritmo visto nos slides em aula.
        if (codigoP0 == 0 && codigoP1 == 0)
            return 1;

        if (verificaBit(codigoP0, codigoP1))
            return 0;

        if (codigoP0 != 0) {
            if (codigoP0 & 1) {
                reta.inicio.x = xmin;
                reta.inicio.y += (xmin - reta.inicio.x) * (reta.fim.y - reta.inicio.y) / (reta.fim.x - reta.inicio.x);
            } else if (codigoP0 & 2) {
                reta.inicio.x = xmax;
                reta.inicio.y += (xmax - reta.inicio.x) * (reta.fim.y - reta.inicio.y) / (reta.fim.x - reta.inicio.x);
            } else if (codigoP0 & 4) {
                reta.inicio.y = ymin;
                reta.inicio.x += (ymin - reta.inicio.y) * (reta.fim.x - reta.inicio.x) / (reta.fim.y - reta.inicio.y);
            } else if (codigoP0 & 8) {
                reta.inicio.y = ymax;
                reta.inicio.x += (ymax - reta.inicio.y) * (reta.fim.x - reta.inicio.x) / (reta.fim.y - reta.inicio.y);
            }
            // Atualiza código do ponto inicial.
            codigoP0 = codificarPonto(reta.inicio, xmin, xmax, ymin, ymax);
        } else if (codigoP1 != 0) {
            if (codigoP1 & 1) {
                reta.fim.x = xmin;
                reta.fim.y += (xmin - reta.fim.x) * (reta.inicio.y - reta.fim.y) / (reta.inicio.x - reta.fim.x);
            } else if (codigoP1 & 2) {
                reta.fim.x = xmax;
                reta.fim.y += (xmax - reta.fim.x) * (reta.inicio.y - reta.fim.y) / (reta.inicio.x - reta.fim.x);
            } else if (codigoP1 & 4) {
                reta.fim.y = ymin;
                reta.fim.x += (ymin - reta.fim.y) * (reta.inicio.x - reta.fim.x) / (reta.inicio.y - reta.fim.y);
            } else if (codigoP1 & 8) {
                reta.fim.y = ymax;
                reta.fim.x += (ymax - reta.fim.y) * (reta.inicio.x - reta.fim.x) / (reta.inicio.y - reta.fim.y);
            }
            // Atualizando código do ponto final.
            codigoP1 = codificarPonto(reta.fim, xmin, xmax, ymin, ymax);
        }
    }
}
// Função que percorre um conjunto de retas e retorna o índice da primeira reta dentro da tolerância t.
int retornaReta(float mx, float my, int tolerancia) {
    for (int i = 0; i < quantidade_retas; i++) {
        if (selecionarReta(retas[i],tolerancia,mx,my)){
            return i;// Retorna o índice da primeira reta dentro da tolerância.
        }
    }
    return -1;// Nenhuma reta encontrada dentro da tolerância.
}
// Função que verifica uma reta dentro da tolerância em relação à posição (mx, my) e a remove se encontrada.
int retaSelecionada(float mx, float my, int tolerancia) {
    for (int i = 0; i < quantidade_retas; i++) {
        if (selecionarReta(retas[i], tolerancia, mx, my)) { // Remove a reta encontrada e ajusta vetor.
            for (int j = i; j < quantidade_retas - 1; j++) {
                retas[j] = retas[j + 1];
            }
            quantidade_retas--; // Atualiza a quantidade de retas após a remoção.

            return 1; // Retorna 1 indicando que uma reta foi removida.
        }
    }

    return 0; // Retorna 0 indicando que nenhuma reta foi removida.
}
// Função que verifica se um ponto (mx, my) selecionado está dentro de algum polígono com base no algoritmo visto nos slides em aula.
int poligonoSelecionado(float mx, float my) {
    for (int i = 0; i < quantidade_poligonos; i++) {
        Poligono p = poligonos[i];
        int intersecoes = 0;

        // Loop para contar interseções entre as arestas do polígono e a linha horizontal no ponto (mx, my).
        for (int j = 0; j < p.quantidade_vertices; j++) {
            Ponto p1 = p.vertices[j];
            Ponto p2 = p.vertices[(j + 1) % p.quantidade_vertices];

            if ((p1.y > my) && (p2.y > my)) continue;
            if ((p1.y < my) && (p2.y < my)) continue;
            if ((p1.x < mx) && (p2.x < mx)) continue;

            // Arestas à direita, uma acima e outra abaixo do ponto em y.
            if (((p1.x > mx) && (p2.x > mx)) && (((p1.y > my) && (p2.y < my)) || ((p1.y < my) && (p2.y > my)))) {
                intersecoes++;
                continue;
            }

            // Caso especial: Aresta horizontal em cima do ponto.
            if (p1.y == p2.y) continue;

            // Calcula a abscissa da interseção.
            float xi = p1.x + (my - p1.y) * (p2.x - p1.x) / (p2.y - p1.y);

            // Verificação com interseção.
            if (xi > mx) {
                intersecoes++;
            }
        }
        // Verifica se o número de interseções é ímpar, indicando que o ponto está dentro do polígono.
        if (intersecoes % 2 == 1) return i;// Polígono selecionado
    }


    return -1;// Nenhum polígono selecionado
}
// Função que verifica um polígono selecionado e o remove atualizando o vetor.
int poligonoExcluido(float mx, float my) {
    int i = poligonoSelecionado(mx,my);
    if(i == -1) return 0; // Nenhum polígono selecionado, nenhum polígono removido.
    for (int j = i; j < quantidade_poligonos - 1; j++) { // Remove o polígono encontrado e atualiza vetor.
            poligonos[j] = poligonos[j + 1];
    }
    quantidade_poligonos--;// Atualiza a quantidade de polígonos após a remoção.
    return 1;// Retorna 1 indicando que um polígono foi removido.
}

// Função para aplicar a translação a um ponto baseado nos slides vistos em aula.
Ponto transladar_ponto(Ponto ponto, float tx, float ty) {
    float matrizTranslacao[3][3] = {
        {1, 0, tx},
        {0, 1, ty},
        {0, 0, 1}
    };

    float matrizPonto[3][1] = {
        {ponto.x},
        {ponto.y},
        {1}
    };

    float resultado[3][1];
    multiplicarMatrizes(resultado, matrizTranslacao, matrizPonto);

    ponto.x = resultado[0][0];
    ponto.y = resultado[1][0];

    return ponto;
}

// Função para aplicar a translação a uma reta.
Reta transladar_reta(Reta reta, float tx, float ty) {
    reta.inicio = transladar_ponto(reta.inicio, tx, ty);
    reta.fim = transladar_ponto(reta.fim, tx, ty);
    return reta;
}

// Função para aplicar a translação a um polígono.
Poligono transladar_poligono(Poligono poligono, float tx, float ty) {
    for (int i = 0; i < poligono.quantidade_vertices; i++) {
        poligono.vertices[i] = transladar_ponto(poligono.vertices[i], tx, ty);
    }
    return poligono;
}
// Função para aplicar a rotação a um ponto baseado nos slides vistos em aula.
Ponto rotacionar_ponto(Ponto ponto, float teta) {
    float matrizRotacao[3][3] = {
        {cos(teta), -sin(teta), 0},
        {sin(teta), cos(teta), 0},
        {0, 0, 1}
    };

    float matrizPonto[3][1] = {
        {ponto.x},
        {ponto.y},
        {1}
    };

    float resultado[3][1];
    multiplicarMatrizes(resultado, matrizRotacao, matrizPonto);

    ponto.x = resultado[0][0];
    ponto.y = resultado[1][0];

    return ponto;
}
// Função para aplicar a rotação a uma reta.
Reta rotacionar_reta(Reta reta, float mx, float my, Ponto pontoinic){
    Ponto centroide;
    centroide.x = (reta.fim.x + reta.inicio.x) / 2;
    centroide.y = (reta.fim.y + reta.inicio.y) / 2;
    float teta = atan2(my - centroide.y,mx - centroide.x) - atan2(pontoinic.y - centroide.y,pontoinic.x - centroide.x);// Calcula o ângulo de rotação em relação ao ponto (mx, my).
    reta.inicio = transladar_ponto(reta.inicio, -centroide.x,-centroide.y);
    reta.fim = transladar_ponto(reta.fim, -centroide.x,-centroide.y);
    reta.inicio = rotacionar_ponto(reta.inicio, teta);
    reta.fim = rotacionar_ponto(reta.fim, teta);
    reta.inicio = transladar_ponto(reta.inicio, centroide.x,centroide.y);
    reta.fim = transladar_ponto(reta.fim, centroide.x,centroide.y);
    return reta;
}
// Função para aplicar a rotação a um polígono.
Poligono rotacionar_poligono(Poligono poligono, float mx, float my, Ponto pontoinic){
    Ponto centroide;
    centroide.x = 0;
    centroide.y = 0;
    float teta;
    for(int i = 0; i < poligono.quantidade_vertices;i++){
        centroide.x += poligono.vertices[i].x;
        centroide.y += poligono.vertices[i].y;
    }
    centroide.x = centroide.x / poligono.quantidade_vertices;
    centroide.y = centroide.y / poligono.quantidade_vertices;
    teta = atan2(my - centroide.y,mx - centroide.x) - atan2(pontoinic.y - centroide.y,pontoinic.x - centroide.x);// Calcula o ângulo de rotação em relação ao ponto (mx, my).
    for(int i = 0; i < poligono.quantidade_vertices; i++){
        poligono.vertices[i] = transladar_ponto(poligono.vertices[i], -centroide.x,-centroide.y);
        poligono.vertices[i] = rotacionar_ponto(poligono.vertices[i], teta);
        poligono.vertices[i] = transladar_ponto(poligono.vertices[i], centroide.x,centroide.y);
    }
    return poligono;
}
// Função baseada nos slides vistos em sala, para aplicar o escalonamento nas retas e poligonos
Ponto escalonar(Ponto ponto, float sx, float sy) {
    float matrizEscala[3][3] = {
        {sx, 0, 0},
        {0, sy, 0},
        {0, 0, 1}
    };

    float matrizPonto[3][1] = {
        {ponto.x},
        {ponto.y},
        {1}
    };

    float resultado[3][1];
    multiplicarMatrizes(resultado, matrizEscala, matrizPonto);

    ponto.x = resultado[0][0];
    ponto.y = resultado[1][0];

    return ponto;
}
// Função para aplicar o escalonamento a uma reta.
Reta escalonar_reta(Reta reta, float mx, float my, float pontoinic_x, float pontoinic_y){
    Ponto centroide;// Calcula o centroide.
    centroide.x = (reta.fim.x + reta.inicio.x) / 2;
    centroide.y = (reta.fim.y + reta.inicio.y) / 2;
    // Calcula as projeções do mouse e do ponto inicial na reta, fiz isso baseado em produto escalar e vetores para basicamente fazer a reta se limitar a
    // escalar apenas na direção dela própria, pois sem isso, devido a ser controlado pelo mouse, a reta estava se movimentando enquanto escalava.
    float delta_x = reta.fim.x - reta.inicio.x;
    float delta_y = reta.fim.y - reta.inicio.y;
    float projecao_mouse = ((((mx - centroide.x) * delta_x) + ((my - centroide.y) * delta_y)) / (delta_x*delta_x+delta_y*delta_y));
    float projecao_inicial = ((((pontoinic_x - centroide.x) * delta_x) + ((pontoinic_y - centroide.y) * delta_y))  / (delta_x*delta_x+delta_y*delta_y));
    float sx = abs(projecao_mouse / projecao_inicial);
    float sy = abs(projecao_mouse / projecao_inicial);
    if(pontoinic_x == centroide.x || mx == centroide.x ){
        sx = 1;
    }
    if(pontoinic_y == centroide.y || my == centroide.y){
        sy = 1;
    }
    float s = sqrt((((sx*sx) + (sy*sy)))/2);

    reta.inicio = transladar_ponto(reta.inicio, -centroide.x,-centroide.y);
    reta.fim = transladar_ponto(reta.fim, -centroide.x,-centroide.y);
    reta.inicio = escalonar(reta.inicio, s, s);
    reta.fim = escalonar(reta.fim, s, s);
    reta.inicio = transladar_ponto(reta.inicio, centroide.x,centroide.y);
    reta.fim = transladar_ponto(reta.fim, centroide.x,centroide.y);
    return reta;
}
// Função para aplicar o escalonamento a um polígono.
Poligono escalonar_poligono(Poligono poligono, float mx, float my, float pontoinic_x, float pontoinic_y){
    Ponto centroide;
    centroide.x = 0;
    centroide.y = 0;
    for(int i = 0; i < poligono.quantidade_vertices;i++){
        centroide.x += poligono.vertices[i].x;
        centroide.y += poligono.vertices[i].y;
    }
    centroide.x = centroide.x / poligono.quantidade_vertices;
    centroide.y = centroide.y / poligono.quantidade_vertices;
    float sx = abs((mx - centroide.x) / (pontoinic_x - centroide.x));
    float sy = abs((my - centroide.y) / (pontoinic_y - centroide.y));
    for (int i = 0; i < poligono.quantidade_vertices; i++) {
        poligono.vertices[i] = transladar_ponto(poligono.vertices[i], -centroide.x, -centroide.y);
        poligono.vertices[i] = escalonar(poligono.vertices[i], sx, sy);
        poligono.vertices[i] = transladar_ponto(poligono.vertices[i], centroide.x, centroide.y);
    }

    return poligono;
}

// Função para aplicar a transformação de reflexão baseado nos slides vistos em sala.
Ponto refletir(Ponto ponto, int eixo) {
    float matrizPonto[3][1] = {{ponto.x},{ponto.y},{1}};
    float matrizReflexao_em_y[3][3] = {{1, 0, 0},{0, -1, 0},{0, 0, 1}};
    float matrizReflexao_em_x[3][3] = {{-1, 0, 0},{0, 1, 0},{0, 0, 1}};
    float matrizReflexaoXY[3][3] = {{-1, 0, 0},{0, -1, 0},{0, 0, 1}};
    float matriz[3][3];
    if(eixo == 0) copiarMatrizes(matriz,matrizReflexao_em_y);
    else if(eixo == 1) copiarMatrizes(matriz,matrizReflexao_em_x);
    else if(eixo == 2) copiarMatrizes(matriz,matrizReflexaoXY);
    float resultado[3][1];
    multiplicarMatrizes(resultado, matriz, matrizPonto);
    ponto.x = resultado[0][0];
    ponto.y = resultado[1][0];
    return ponto;
}
// Função para aplicar a reflexão a uma reta.
Reta refletir_reta(Reta reta, int eixo){
    Ponto centroide;
    centroide.x = (reta.fim.x + reta.inicio.x) / 2;
    centroide.y = (reta.fim.y + reta.inicio.y) / 2;
    //reflete em relação ao centro do objeto se "descomentar".
    //reta.inicio = transladar_ponto(reta.inicio, -centroide.x,-centroide.y);
    //reta.fim = transladar_ponto(reta.fim, -centroide.x,-centroide.y);
    reta.inicio = refletir(reta.inicio, eixo);
    reta.fim = refletir(reta.fim, eixo);
    //reta.inicio = transladar_ponto(reta.inicio, centroide.x,centroide.y);
    //reta.fim = transladar_ponto(reta.fim, centroide.x,centroide.y);
    return reta;
}
// Função para aplicar a reflexão a um polígono.
Poligono refletir_poligono(Poligono poligono, int eixo){
    Ponto centroide;
    centroide.x = 0;
    centroide.y = 0;
    for(int i = 0; i < poligono.quantidade_vertices;i++){
        centroide.x += poligono.vertices[i].x;
        centroide.y += poligono.vertices[i].y;
    }
    centroide.x = centroide.x / poligono.quantidade_vertices;
    centroide.y = centroide.y / poligono.quantidade_vertices;
    for (int i = 0; i < poligono.quantidade_vertices; i++) {
        //reflete em relação ao centro do objeto se "descomentar".
        //poligono.vertices[i] = transladar_ponto(poligono.vertices[i], -centroide.x, -centroide.y);
        poligono.vertices[i] = refletir(poligono.vertices[i], eixo);
        //poligono.vertices[i] = transladar_ponto(poligono.vertices[i], centroide.x, centroide.y);
    }

    return poligono;
}

// Função para aplicar a transformação de cisalhamento baseado nos slides vistos em sala.
Ponto cisalhar(Ponto ponto, float sh, int eixo) {
    float matrizcisalha_em_x[3][3] = {
        {1, sh, 0},
        {0, 1, 0},
        {0, 0, 1}
    };
    float matrizcisalha_em_y[3][3] = {
        {1, 0, 0},
        {sh, 1, 0},
        {0, 0, 1}
    };

    float matrizPonto[3][1] = {
        {ponto.x},
        {ponto.y},
        {1}
    };

    if(eixo == 0){
        float resultado[3][1];
        multiplicarMatrizes(resultado, matrizcisalha_em_x, matrizPonto);
        ponto.x = resultado[0][0];
        ponto.y = resultado[1][0];
    }
    if(eixo == 1){
        float resultado[3][1];
        multiplicarMatrizes(resultado, matrizcisalha_em_y, matrizPonto);
        ponto.x = resultado[0][0];
        ponto.y = resultado[1][0];
    }

    return ponto;
}
// Função para aplicar cisalhamento a uma reta.
Reta cisalhar_reta(Reta reta, float mx, float my, Ponto pontoinic, int eixo){
    Ponto centroide;
    centroide.x = (reta.fim.x + reta.inicio.x) / 2;
    centroide.y = (reta.fim.y + reta.inicio.y) / 2;

    float shx = (mx - pontoinic.x) / (pontoinic.y - centroide.y);
    float shy = (my - pontoinic.y) / (pontoinic.x - centroide.x);
    if(my == centroide.y) shx = 0;
    if(mx == centroide.x) shy = 0;
    if(eixo == 0){
        reta.inicio = transladar_ponto(reta.inicio, -centroide.x,-centroide.y);
        reta.fim = transladar_ponto(reta.fim, -centroide.x,-centroide.y);
        reta.inicio = cisalhar(reta.inicio, shx ,eixo);
        reta.fim = cisalhar(reta.fim, shx ,eixo);
        reta.inicio = transladar_ponto(reta.inicio, centroide.x,centroide.y);
        reta.fim = transladar_ponto(reta.fim, centroide.x,centroide.y);
    }
    if(eixo == 1){
        reta.inicio = transladar_ponto(reta.inicio, -centroide.x,-centroide.y);
        reta.fim = transladar_ponto(reta.fim, -centroide.x,-centroide.y);
        reta.inicio = cisalhar(reta.inicio, shy ,eixo);
        reta.fim = cisalhar(reta.fim, shy ,eixo);
        reta.inicio = transladar_ponto(reta.inicio, centroide.x,centroide.y);
        reta.fim = transladar_ponto(reta.fim, centroide.x,centroide.y);
    }
    return reta;
}
// Função para aplicar cisalhamento a um polígono.
Poligono cisalhar_poligono(Poligono poligono, float mx, float my, Ponto pontoinic, int eixo){
    Ponto centroide;
    centroide.x = 0;
    centroide.y = 0;
    for(int i = 0; i < poligono.quantidade_vertices;i++){
        centroide.x += poligono.vertices[i].x;
        centroide.y += poligono.vertices[i].y;
    }
    centroide.x = centroide.x / poligono.quantidade_vertices;
    centroide.y = centroide.y / poligono.quantidade_vertices;

    float shx = (mx - pontoinic.x) / (pontoinic.y - centroide.y);
    float shy = (my - pontoinic.y) / (pontoinic.x - centroide.x);
    if(my == centroide.y) shx = 0;
    if(mx == centroide.x) shy = 0;

    if(eixo == 0){
        for(int i = 0; i < poligono.quantidade_vertices; i++){
            poligono.vertices[i] = transladar_ponto(poligono.vertices[i], -centroide.x,-centroide.y);
            poligono.vertices[i] = cisalhar(poligono.vertices[i], shx, eixo);
            poligono.vertices[i] = transladar_ponto(poligono.vertices[i], centroide.x,centroide.y);
        }
    }
    if(eixo == 1){
        for(int i = 0; i < poligono.quantidade_vertices; i++){
            poligono.vertices[i] = transladar_ponto(poligono.vertices[i], -centroide.x,-centroide.y);
            poligono.vertices[i] = cisalhar(poligono.vertices[i], shy, eixo);
            poligono.vertices[i] = transladar_ponto(poligono.vertices[i], centroide.x,centroide.y);
        }
    }

    return poligono;
}
// Função salva as informações da tela em um arquivo chamado "arquivo_salvo" onde arquivo contém basicamente
// a quantidade de pontos, retas e polígonos e seus dados das estruturas.
void salvar_arquivo(){
    ofstream arq("arquivo_salvo", ios::trunc|ios::out);
    if(!arq){// Verifica se o arquivo foi criado corretamente.
        cout << "Arquivo não criado!";
        return;
    }else{
        cout << "Arquivo criado!";
    }// Escreve a quantidade de pontos, retas e polígonos no arquivo.
    arq << quantidade_pontos <<" "<< quantidade_retas <<" "<< quantidade_poligonos << "\n";
    // Escreve os dados dos pontos no arquivo.
    for(int i = 0; i < quantidade_pontos;i++){
        arq << pontos[i].x << " " << pontos[i].y << " ";
        for(int c = 0; c < 3; c++){
            arq << pontos[i].cor[c] << " ";
        }
    }
    arq << "\n";// Escreve os dados das retas no arquivo.
    for(int j = 0; j < quantidade_retas;j++){
            arq << retas[j].inicio.x << " " << retas[j].inicio.y << " ";
            for(int c = 0; c < 3; c++){
                arq << retas[j].inicio.cor[c] << " ";
            }
            arq << retas[j].fim.x << " " << retas[j].fim.y << " ";
            for(int c = 0; c < 3; c++){
                arq << retas[j].fim.cor[c] << " ";
            }
    }
    arq << "\n";// Escreve os dados dos poligonos no arquivo.
    for(int k = 0; k < quantidade_poligonos;k++){
        arq << poligonos[k].quantidade_vertices << " ";
        for(int c = 0; c < 3; c++){
                arq << poligonos[k].cor[c] << " ";
            }
        for (int t = 0; t < poligonos[k].quantidade_vertices; t++){
            arq << poligonos[k].vertices[t].x << " " << poligonos[k].vertices[t].y << " ";
            for(int c = 0; c < 3; c++){
                arq << poligonos[k].vertices[t].cor[c] << " ";
            }
        }
    }
    arq << "\n";
    arq.close();
    return;
}
// Função carrega as informações do arquivo chamado "arquivo_salvo" para o programa.
void carregar_arquivo(){
    ifstream arq("arquivo_salvo",ios::in);
    if(!arq){// Verifica se o arquivo foi carregado corretamente.
        cout << "Arquivo não carregado!";
        return;
    }else{
        cout << "Arquivo carregado!";
    } // Lê a quantidade de pontos, retas e polígonos do arquivo.
    arq >> quantidade_pontos >> quantidade_retas >> quantidade_poligonos ;
    for(int i = 0; i < quantidade_pontos;i++){// Lê os dados dos pontos do arquivo.
        arq >> pontos[i].x >> pontos[i].y;
        for(int c = 0; c < 3; c++){
            arq >> pontos[i].cor[c];
        }
    }
    for(int j = 0; j < quantidade_retas;j++){ // Lê os dados das retas do arquivo.
            arq >> retas[j].inicio.x >> retas[j].inicio.y;
            for(int c = 0; c < 3; c++){
                arq >> retas[j].inicio.cor[c];
            }
            arq >> retas[j].fim.x >> retas[j].fim.y;
            for(int c = 0; c < 3; c++){
                arq >> retas[j].fim.cor[c];
            }
    }
    for(int k = 0; k < quantidade_poligonos;k++){ // Lê os dados dos poligonos do arquivo.
        arq >> poligonos[k].quantidade_vertices;
        for(int c = 0; c < 3; c++){
                arq >> poligonos[k].cor[c];
            }
        for (int t = 0; t < poligonos[k].quantidade_vertices; t++){
            arq >> poligonos[k].vertices[t].x >> poligonos[k].vertices[t].y;
            for(int c = 0; c < 3; c++){
                arq >> poligonos[k].vertices[t].cor[c];
            }
        }
    }
    arq.close();
    glutPostRedisplay();
    return;
}
// Função display padrão.
void display() {
    glClear(GL_COLOR_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    desenharPontos();
    desenharRetas();
    desenharPoligonos();
    glFlush();
}
// Função de todos os comandos do teclado.
void retorno_teclado(unsigned char key, int x, int y) {
    switch (key) {
        case 27: // Código ASCII para a tecla 'ESC'
            exit(0);
            break;
        case 'p':
            modo_desenho = 0; // Ponto
            break;
        case 'l':
            modo_desenho = 1; // Reta
            break;
        case 'm':
            modo_desenho = 2; // Poligono
            break;
        case 'n': // Completa e Preenche Poligono
            if (modo_desenho == 2) {
                if (desenhaPoligono) {
                    desenhaPoligono = 0;
                    for(int i = 0; i < quantidadePontosPoligono; i++){
                        pontoSelecionado(pontosPoligono[i].x,pontosPoligono[i].y,0.05);
                    }
                    adicionarPoligono(pontosPoligono, quantidadePontosPoligono, poligonos[0].cor);
                    quantidadePontosPoligono = 0;
                    glutPostRedisplay();
                } else {
                    desenhaPoligono = 1;
                }
            }
            break;
        case 't': // translada
            transformacao = 1;
            break;
        case 'r': // rotacao
            transformacao = 2;
            break;
        case 'e': // escala
            transformacao = 3;
            break;
        case 'y': // reflexao y
            transformacao = 4;
            break;
        case 'x': // reflexao x
            transformacao = 5;
            break;
        case 'z': // reflexao xy
            transformacao = 6;
            break;
        case 'c': // cisalhar x
            transformacao = 7;
            break;
        case 'v': // cisalhar y
            transformacao = 8;
        case 's':
            salvar_arquivo();
            break;
        case 'd':
            carregar_arquivo();
            break;
    }
}
// Função de todos os comandos de clique do mouse.
void retorno_mouse(int button, int state, int x, int y) {
    static Ponto pontoInicialLinha;
    static int clicouPrimeiroPonto = 0;
    float mx = (x / 800.0) * 2 - 1;// Tive problemas com a inversão comentada em aula, e não estava conseguindo movimentar o mouse corretamente,
    float my = 1 - (y / 800.0) * 2;// até pesquisar e aplicar essa fórmula e funcionar para o meu código, adptando para o tamanho do meu programa.
    if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN) {
        if (modo_desenho == 0) { // Ponto
            if (transformacao != 0) {
                objetoatual = retornaPonto(mx, my, 0.05);
                pontoInicialTransformacao.x = mx;
                pontoInicialTransformacao.y = my;
            }else {
                float corPonto[3] = {1, 1, 0};
                adicionarPonto(mx, my, corPonto);
            }
        } else if (modo_desenho == 1) { // Reta
            if (transformacao != 0) {
                objetoatual = retornaReta(mx, my, 0.05);
                pontoInicialTransformacao.x = mx;
                pontoInicialTransformacao.y = my;
                if(objetoatual != -1){
                    if(transformacao == 4){ // Reflexão Y
                        retas[objetoatual] = refletir_reta(retas[objetoatual], 0);
                    }
                    if(transformacao == 5){ // Reflexão X
                        retas[objetoatual] = refletir_reta(retas[objetoatual], 1);
                    }
                    if(transformacao == 6){ // Reflexão XY
                        retas[objetoatual] = refletir_reta(retas[objetoatual], 2);
                    }
                }
            }else {
                if (!clicouPrimeiroPonto) {
                    float corReta[3] = {0, 1, 0};
                    adicionarPonto(mx, my, corReta);
                    pontoInicialLinha.x = mx;
                    pontoInicialLinha.y = my;
                    clicouPrimeiroPonto = 1;
                } else {
                    float corReta[3] = {0, 1, 0};
                    adicionarReta(pontoInicialLinha.x, pontoInicialLinha.y, mx, my, corReta);
                    clicouPrimeiroPonto = 0;
                    pontoSelecionado(pontoInicialLinha.x, pontoInicialLinha.y, 0.05);
                }
            }
        } else if (modo_desenho == 2) { // Poligono
                if (transformacao != 0) {
                    objetoatual = poligonoSelecionado(mx, my);
                    pontoInicialTransformacao.x = mx;
                    pontoInicialTransformacao.y = my;
                    if(objetoatual != -1){
                        if(transformacao == 4){ // Reflexão Y
                            poligonos[objetoatual] = refletir_poligono(poligonos[objetoatual], 0);
                        }
                        if(transformacao == 5){ // Reflexão X
                            poligonos[objetoatual] = refletir_poligono(poligonos[objetoatual], 1);
                        }
                        if(transformacao == 6){ // Reflexão XY
                            poligonos[objetoatual] = refletir_poligono(poligonos[objetoatual], 2);
                        }
                    }
                } else {
                    desenhaPoligono = 1;
                    float corPontoPoligono[3] = {1, 0, 0};
                    adicionarPonto(mx, my, corPontoPoligono);
                    pontosPoligono[quantidadePontosPoligono] = pontos[quantidade_pontos - 1];
                    quantidadePontosPoligono++;
                }
            }
    }else if (button == GLUT_LEFT_BUTTON && state == GLUT_UP) { // Solta clique
        transformacao = 0;
    } else if (button == GLUT_RIGHT_BUTTON && state == GLUT_DOWN) {
        switch (modo_desenho) {
            case 0:// exclui o ponto
                if (pontoSelecionado(mx, my, 0.05)) {
                    glutPostRedisplay();
                }
                break;
            case 1:// exclui a reta
                if (retaSelecionada(mx, my, 0.05)) {
                    glutPostRedisplay();
                }
                break;
            case 2: // exclui o Polígono
                if (poligonoExcluido(mx, my)) {
                    glutPostRedisplay();
                }
                break;
        }
    }
    glutPostRedisplay();
}
// Função de todos os comandos de "arrasto" ou movimento do mouse.
void retorno_movimento_mouse(int x, int y) {
    float mx = (x / 800.0) * 2 - 1;
    float my = 1 - (y / 800.0) * 2;

    if (transformacao == 1) { // translação
        if (modo_desenho == 0) {
            pontos[objetoatual] = transladar_ponto(pontos[objetoatual], mx - pontoInicialTransformacao.x, my - pontoInicialTransformacao.y);
            pontoInicialTransformacao.x = mx;
            pontoInicialTransformacao.y = my;
        } else if (modo_desenho == 1) {
            if(objetoatual != -1){
                retas[objetoatual] = transladar_reta(retas[objetoatual], mx - pontoInicialTransformacao.x, my - pontoInicialTransformacao.y);
                pontoInicialTransformacao.x = mx;
                pontoInicialTransformacao.y = my;
            }
        } else if (modo_desenho == 2) {
            if(objetoatual != -1){
                poligonos[objetoatual] = transladar_poligono(poligonos[objetoatual], mx - pontoInicialTransformacao.x, my - pontoInicialTransformacao.y);
                pontoInicialTransformacao.x = mx;
                pontoInicialTransformacao.y = my;
            }
        }
    }else if(transformacao == 2){ // Rotação
        if(modo_desenho == 0){
            pontos[objetoatual] = rotacionar_ponto(pontos[objetoatual],atan2(my,mx) - atan2(pontoInicialTransformacao.y,pontoInicialTransformacao.x));
            pontoInicialTransformacao.x = mx;
            pontoInicialTransformacao.y = my;
        }else if(modo_desenho == 1){
            if(objetoatual != -1){
                retas[objetoatual] = rotacionar_reta(retas[objetoatual],mx,my,pontoInicialTransformacao);
                pontoInicialTransformacao.x = mx;
                pontoInicialTransformacao.y = my;
            }
        }else if(modo_desenho == 2){
            if(objetoatual != -1){
                poligonos[objetoatual] = rotacionar_poligono(poligonos[objetoatual],mx,my,pontoInicialTransformacao);
                pontoInicialTransformacao.x = mx;
                pontoInicialTransformacao.y = my;
            }
        }
    }else if (transformacao == 3) { // Escalar
        if (modo_desenho == 1) {
            if(objetoatual != -1){
                retas[objetoatual] = escalonar_reta(retas[objetoatual], mx, my, pontoInicialTransformacao.x,pontoInicialTransformacao.y);
                pontoInicialTransformacao.x = mx;
                pontoInicialTransformacao.y = my;
            }
        } else if (modo_desenho == 2) {
            if(objetoatual != -1){
                poligonos[objetoatual] = escalonar_poligono(poligonos[objetoatual], mx, my, pontoInicialTransformacao.x,pontoInicialTransformacao.y);
                pontoInicialTransformacao.x = mx;
                pontoInicialTransformacao.y = my;
            }
        }
    }else if (transformacao == 7) { // Cisalhar X
        if (modo_desenho == 1){
            if(objetoatual != -1){
                retas[objetoatual] = cisalhar_reta(retas[objetoatual], mx, my, pontoInicialTransformacao,0);
                pontoInicialTransformacao.x = mx;
                pontoInicialTransformacao.y = my;
            }
        }else if (modo_desenho == 2){
            if(objetoatual != -1){
                poligonos[objetoatual] = cisalhar_poligono(poligonos[objetoatual], mx, my, pontoInicialTransformacao,0);
                pontoInicialTransformacao.x = mx;
                pontoInicialTransformacao.y = my;
            }
        }
    }else if (transformacao == 8) { // Cisalhar Y
        if (modo_desenho == 1){
            if(objetoatual != -1){
                retas[objetoatual] = cisalhar_reta(retas[objetoatual], mx, my, pontoInicialTransformacao,1);
                pontoInicialTransformacao.x = mx;
                pontoInicialTransformacao.y = my;
            }
        }else if (modo_desenho == 2){
            if(objetoatual != -1){
                poligonos[objetoatual] = cisalhar_poligono(poligonos[objetoatual], mx, my, pontoInicialTransformacao,1);
                pontoInicialTransformacao.x = mx;
                pontoInicialTransformacao.y = my;
            }
        }
    }
    glutPostRedisplay();
}
// Função init padrão.
void init() {
    glClearColor(0.0, 0.15, 0.25, 1.0); // Define a cor de fundo
    glMatrixMode(GL_PROJECTION); // Carrega a matriz de projeção
    glOrtho(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0); // Define a projeção ortogonal 2D
    glutKeyboardFunc(retorno_teclado);
}
int main(int argc, char** argv) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
    glutInitWindowSize(800, 800);
    glutInitWindowPosition(100, 100);
    glutCreateWindow("Pontos, Retas e Poligonos com OpenGL");
    init();
    glutMotionFunc(retorno_movimento_mouse);
    glutMouseFunc(retorno_mouse);
    glutKeyboardFunc(retorno_teclado);
    glutDisplayFunc(display);
    glutMainLoop();
    return 0;
}
