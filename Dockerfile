FROM ubuntu:22.04

# Evitar interacciones al instalar paquetes
ENV DEBIAN_FRONTEND=noninteractive

# Instalar dependencias de compilación y ejecución
RUN apt-get update && apt-get install -y \
    git \
    build-essential \
    g++ \
    make \
    cmake \
    unzip \
    python3 \
    python3-pip \
    libgomp1 \
    libgfortran5 \
    perl \
    && rm -rf /var/lib/apt/lists/*

# Crear directorios necesarios
RUN mkdir -p /opt/unidesign \
    && mkdir -p /workspace

# Directorio de trabajo principal
WORKDIR /opt/unidesign

# Clonar el repositorio UniDesign
RUN git clone "https://github.com/XRed8X/UniDesign.git" source

# Compilar UniDesign
WORKDIR /opt/unidesign/source
RUN chmod +x build.sh && ./build.sh

# IMPORTANTE: Dar permisos de ejecución a TODOS los archivos ejecutables
RUN chmod -R +x /opt/unidesign/source/evolution/bin/ || true && \
    chmod -R +x /opt/unidesign/source/extbin/ || true && \
    chmod +x /opt/unidesign/source/*.pl || true && \
    chmod +x /opt/unidesign/source/*.sh || true && \
    chmod +x /opt/unidesign/source/wread/*.pl || true

# Verificar que el binario principal existe y darle permisos
RUN if [ -f /opt/unidesign/source/UniDesign ]; then \
        chmod +x /opt/unidesign/source/UniDesign && \
        ln -s /opt/unidesign/source/UniDesign /usr/local/bin/UniDesign; \
    elif [ -f /opt/unidesign/source/bin/UniDesign ]; then \
        chmod +x /opt/unidesign/source/bin/UniDesign && \
        ln -s /opt/unidesign/source/bin/UniDesign /usr/local/bin/UniDesign; \
    else \
        echo "ERROR: UniDesign binary not found after build!" && exit 1; \
    fi

# Añadir directorios al PATH
ENV PATH="/opt/unidesign/source:/opt/unidesign/source/bin:/opt/unidesign/source/extbin:/opt/unidesign/source/evolution/bin:${PATH}"

# Directorio de trabajo para el usuario
WORKDIR /workspace

# Crear entrypoint script
RUN echo '#!/bin/bash' > /opt/unidesign/entrypoint.sh && \
    echo 'exec /opt/unidesign/source/UniDesign "$@"' >> /opt/unidesign/entrypoint.sh && \
    chmod +x /opt/unidesign/entrypoint.sh

# Definir entrypoint
ENTRYPOINT ["/opt/unidesign/entrypoint.sh"]